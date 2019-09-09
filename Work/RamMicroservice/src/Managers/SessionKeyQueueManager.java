package Managers;

import QueueClass.ResQueueItem;
import QueueClass.SessionArrValue;
import QueueClass.SessionArrValueConverter;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Timestamp;
import java.util.Date;
import java.util.Map;

public class SessionKeyQueueManager {
    QueueManager queue;
    LogManager log_manager;
    Gson gson;
    
    String name_class = "SessionKeyQueueManager";
    
    public SessionKeyQueueManager(LogManager log_manager) {
        this.log_manager = log_manager;
        queue = new QueueManager(log_manager, "SessionKeyQueue");
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(SessionArrValue.class, new SessionArrValueConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        gson = builder.create();
    }
    
    public ResData addSessionKeyItem(ResQueueItem item) {
        String SessionKey = "";
        try {
            for (Object buffitem : queue.getQueue().values()) {
                if (((SessionArrValue) buffitem).getId() == item.getId()) {
                    queue.DeleteItemFromItem(buffitem);
                    break;
                }
            }

            SessionKey = queue.RandMd5Key();
            queue.AddItem(SessionKey, new SessionArrValue(item.getId(), new Timestamp(System.currentTimeMillis())));
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "addSessionKeyItem(). Exception", name_class);
        }
        return new ResData(true, gson.toJson(new SQLResult(gson.toJson(new ResQueueItem(item.getId(), SessionKey), ResQueueItem.class), 1, 1), SQLResult.class));
    }
    
    public ResData CheckUserId(String key) {
        if(queue.getQueue().containsKey(key)) {
            SessionArrValue buff = ((SessionArrValue)queue.getItemForKey(key));
            return new ResData(true, gson.toJson(new SQLResult("Ok", buff.getId(), buff.getId())));
        }
        else 
            return new ResData(true, gson.toJson(new SQLResult("not Ok",1,1)));
    }
    
    public ResData GetSessionKeyItems() {
        String list = "";
        try {
            for (Map.Entry<String, Object> item : queue.getQueue().entrySet()) {
                SessionArrValue buff = (SessionArrValue) item.getValue();
                list = list + "Session key: " + item.getKey() + " user_id: " + buff.getId() + " date: " + buff.getDate().toString() + "\n";
            }
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "GetSessionKeyItems(). Exception", name_class);
        }
        return new ResData(true, gson.toJson(new SQLResult(list,1,1), SQLResult.class));
    }
    
    public ResData GetSessionKeyItem(ResQueueItem item) {
        String session_key = "";
        try {
            for (Map.Entry<String, Object> entry : queue.getQueue().entrySet()) {
        	SessionArrValue buff = (SessionArrValue) entry.getValue();
		    if (buff.getId() == item.getId())  {
			session_key = entry.getKey();
			break;
		    }
            }
	    if (session_key.equals("")) {
		return this.addSessionKeyItem(item);
		//new ResData(true, gson.toJson(new SQLResult(gson.toJson(new ResQueueItem(item.getId(), SessionKey), ResQueueItem.class), 1, 1), SQLResult.class));
	    }
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "GetSessionKeyItem(). Exception", name_class);
        }
        return new ResData(true, gson.toJson(new SQLResult(gson.toJson(new ResQueueItem(item.getId(), session_key), ResQueueItem.class), 1, 1), SQLResult.class));
    }
    
    public ResData RefresQueue(Date dateNow, int max_sec) {
        try {
            for (int i = 0; i < queue.getQueue().size(); i++) {
                if (dateNow.getTime() - ((SessionArrValue) (queue.getQueue().values().toArray()[i])).getDate().getTime() > max_sec) {
                    queue.getQueue().remove(queue.getQueue().keySet().toArray()[i]);
                    if (queue.getQueue().isEmpty()) {
                        break;
                    }
                    i = i - 1;
                }
            }
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "RefresQueue(). Exception", name_class);
        }
        return new ResData(true, gson.toJson(new SQLResult("Ok",1,1)));
    }
}
