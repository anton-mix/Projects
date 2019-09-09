package Managers;

import QueueClass.ResQueueItem;
import QueueClass.UserArrValue;
import QueueClass.UserArrValueConverter;
import QueueClass.UserAuthModel;
import QueueClass.UserAuthModelConverter;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Timestamp;
import java.util.Date;
import java.util.Map;

public class UserAuthMessageQueueManager {
    QueueManager queue;
    LogManager log_manager;
    Gson gson;
    
    String name_class = "UserAuthMessageQueueManager";
    
    public UserAuthMessageQueueManager(LogManager log_manager) {
        this.log_manager = log_manager;
        queue = new QueueManager(log_manager, "UserAuthMessageQueue");
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(UserAuthModel.class, new UserAuthModelConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        builder.registerTypeAdapter(UserArrValue.class, new UserArrValueConverter());
        gson = builder.create();
    }
    
    public ResData AddUserAuthMessageItem(UserAuthModel item) {
        String Session_key = "";
        try {
            Session_key = queue.RandMd5Key();
            queue.AddItem(Session_key, new UserArrValue(new Timestamp(System.currentTimeMillis()), item));
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddUserAuthMessage(UserAuthModel). Exception", name_class);
        }
        
        return new ResData(true, gson.toJson(new SQLResult(Session_key, 1, 1), SQLResult.class));
    }
    
    public ResData DeleteUserAuthMessageItem(String key) {
        try {
            queue.DeleteItemFromKey(key);
        }
        catch (Exception ex){
            return log_manager.CallError(ex.toString(), "DeleteUseaAuthMessage(String). Exception", name_class);
        }
        
        return new ResData(true, gson.toJson(new SQLResult("Ok",1,1), SQLResult.class));
    }
    
    public ResData DeleteUserAuthMessageItemForEmail(UserAuthModel item) {
        try {
            for(Object item2: queue.getQueue().values()) {
                if(((UserArrValue)item2).getUser().getEmail().equals(item.getEmail())) {
                    queue.DeleteItemFromItem(item2);
                }
            }
            return new ResData(true, gson.toJson(new SQLResult("Ok",1,1), SQLResult.class));
        }
        catch(Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteUserAuthMessageItemForEmail(String). Exception", name_class);
        }
    }
    
    public ResData GetUserAuthMessageItems() {
        String list = "";
        try {
            for (Map.Entry<String, Object> buffitem : queue.getQueue().entrySet()) {
                UserArrValue item = ((UserArrValue) buffitem.getValue());
                list = list + "Session key: " + buffitem.getKey() + " user_id: " + item.getUser().toString() + " date: " + item.getDate().toString() + "\n";
            }
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "GetUserAuthMessageItems(). Exception", name_class);
        }
        return new ResData(true, gson.toJson(new SQLResult(list,1,1), SQLResult.class));
    }
    
    public ResData CheckAuthKey(String key) {
        return new ResData(true, gson.toJson(new SQLResult((queue.getQueue().containsKey(key))?"Ok":"not Ok",1,1)));
    }
    
    public ResData CheckAuthMessageKeyEmail(String email) {
        boolean checkEmail = true;
        try {
            for (Object item : queue.getQueue().values()) {
                if (((UserArrValue) item).getUser().getEmail().equals(email)) {
                    checkEmail = false;
                    break;
                }
            }
            return new ResData(true, gson.toJson(new SQLResult((checkEmail) ? "Ok" : "not Ok", 1, 1)));
        }
        catch(Exception ex) {
            return log_manager.CallError(ex.toString(), "CheckAuthMessageKeyEmail(SessionKey). Exception", name_class);
        }
    }
    
    public ResData CheckAndDeleteUserAuthMessage(String SessionKey) {
        try {
            if (queue.getQueue().containsKey(SessionKey)) {
                UserArrValue buff = ((UserArrValue) queue.getItemForKey(SessionKey));
                DeleteUserAuthMessageItem(SessionKey);
                return new ResData(true, gson.toJson(new SQLResult(gson.toJson(buff.getUser(), UserAuthModel.class) ,1 ,1), SQLResult.class));
            } else {
                return new ResData(true, gson.toJson(new SQLResult(gson.toJson(new UserAuthModel("", "", "", "", "", ""), UserAuthModel.class) ,1 ,1), SQLResult.class));
            }
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "CheckAndDeleteUserAuthMessage(SessionKey). Exception", name_class);
        }
    }
    
    public ResData RefresQueue(Date dateNow, int max_sec) {
        try {
            for (int i = 0; i < queue.getQueue().size(); i++) {
                if (dateNow.getTime() - ((UserArrValue) (queue.getQueue().values().toArray()[i])).getDate().getTime() > max_sec) {
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
        return new ResData(true, gson.toJson(new SQLResult("Ok" ,1 ,1)));
    }
}