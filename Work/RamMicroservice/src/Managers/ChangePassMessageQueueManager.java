package Managers;

import QueueClass.ArrUserPasswordKeyValue;
import QueueClass.ArrUserPasswordKeyValueConverter;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Timestamp;
import java.util.Date;
import java.util.Map;

public class ChangePassMessageQueueManager {
    QueueManager queue;
    LogManager log_manager;
    Gson gson;
    
    String name_class = "ChangePassMessageQueueManager";
    
    public ChangePassMessageQueueManager(LogManager log_manager) {
        this.log_manager = log_manager;
        queue = new QueueManager(log_manager, "ChangePassMessageQueue");
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(ArrUserPasswordKeyValue.class, new ArrUserPasswordKeyValueConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        gson = builder.create();
    }
    
    public ResData AddChangePassMessageItem(ArrUserPasswordKeyValue item) {
        String SessionKey = "";
        try {
            SessionKey = queue.RandMd5Key();
            queue.AddItem(SessionKey, new ArrUserPasswordKeyValue(new Timestamp(System.currentTimeMillis()), item.getLogin()));
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddChangePassMessageItem(ArrUserPasswordKeyValue). Exception", name_class);
        }

        return new ResData(true, gson.toJson(new SQLResult(SessionKey, 1, 1), SQLResult.class));
    }
    
    public ResData DeleteChangePassMessageItem(String key) {
        try {
            queue.DeleteItemFromKey(key);
        }
        catch (Exception ex){
            return log_manager.CallError(ex.toString(), "DeleteChangePassMessageItem(String). Exception", name_class);
        }
        
        return new ResData(true, gson.toJson(new SQLResult("Ok" ,1 ,1), SQLResult.class));
    }
    
    public ResData DeleteChangePassMessageItemFromLogin(ArrUserPasswordKeyValue item) {
        try {
            for(Object item2 : queue.getQueue().values())  {
                if(((ArrUserPasswordKeyValue)item2).getLogin().equals(item.getLogin())) {
                    queue.DeleteItemFromItem(item2);
                    break;
                }
            }
        }
        catch (Exception ex){
            return log_manager.CallError(ex.toString(), "DeleteChangePassMessageItem(String). Exception", name_class);
        }
        
        return new ResData(true, gson.toJson(new SQLResult("Ok",1,1), SQLResult.class));
    }
    
    public ResData GetChangePassMessageItems() {
        String list = "";
        try {
            for (Map.Entry<String, Object> buffitem : queue.getQueue().entrySet()) {
                ArrUserPasswordKeyValue item = ((ArrUserPasswordKeyValue) buffitem.getValue());
                list = list + "Session key: " + buffitem.getKey() + " login: " + item.getLogin() + " date: " + item.getDate().toString() + "\n";
            }
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "GetChangePassMessageItems(String). Exception", name_class);
        }
        return new ResData(true, gson.toJson(new SQLResult(list,1,1), SQLResult.class));
    }
    
    public ResData CheckVhangePassKey(String key) {
        return new ResData(true, gson.toJson(new SQLResult((queue.getQueue().containsKey(key))?"Ok":"not Ok", 1 ,1)));
    }
    
    public ResData CheckChangePassMessageKeyLogin(String login) {
        boolean checkLogin = true;
        for(Object item : queue.getQueue().values()) {
            if(((ArrUserPasswordKeyValue)item).getLogin().equals(login)) {
                checkLogin = false;
                break;
            }
        }
        return new ResData(true, gson.toJson(new SQLResult((checkLogin)?"Ok":"not Ok",1,1)));
    }
    
    public ResData CheckAndDeleteChangePassMessageKey(String SessionKey) {
        try {
            if (queue.getQueue().containsKey(SessionKey)) {
                ArrUserPasswordKeyValue buff = (ArrUserPasswordKeyValue)queue.getItemForKey(SessionKey);
                DeleteChangePassMessageItem(SessionKey);
                return new ResData(true, gson.toJson(new SQLResult(gson.toJson(buff, ArrUserPasswordKeyValue.class),1 ,1), SQLResult.class));
            } else {
                return new ResData(true, gson.toJson(new SQLResult("not Ok" ,1 ,1), SQLResult.class));
            }
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "CheckAndDeleteUserAuthMessage(SessionKey). Exception", name_class);
        }
    }
    
    public ResData RefresQueue(Date dateNow, int max_sec) {
        try {
            for (int i = 0; i < queue.getQueue().size(); i++) {
                if (dateNow.getTime() - ((ArrUserPasswordKeyValue) (queue.getQueue().values().toArray()[i])).getDate().getTime() > max_sec) {
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
