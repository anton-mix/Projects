package services;

import Managers.LogManager;
import Managers.MessageManager;
import MyEnum.RestNotificationLogEnum;
import MyInterface.RestInterface;
import RestObject.RestObject;
import RestObject.RestObjectConverter;
import SQLClass.GlobalVariables;
import SQLClass.ResData;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import javax.ws.rs.Consumes;
import javax.ws.rs.FormParam;
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.core.Response;

@Path("ManagerNotificationLog")
public class ManagerNotificationLogRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerNotificationLogRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerNotificationLogRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerNotificationLogRest");
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            answer = log_manager.CallError(ex.toString(), "ManagerNotificationLogRest(). Init", name_class);
        } 
    }
    
    @POST
    @Path("getNotificationLog")
    @Consumes("application/x-www-form-urlencoded")
    public Response getNotificationLog(@FormParam("json") String json) {
        answer = CreateMessage(RestNotificationLogEnum.Notification_log_select.ordinal(), json, "getNotificationLog(String). Send message", "notification_log_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getCountNotificationLog")
    @Consumes("application/x-www-form-urlencoded")
    public Response getCountNotificationLog(@FormParam("json") String json) {
        answer = CreateMessage(RestNotificationLogEnum.Notification_log_select_count.ordinal(), json, "getCountNotificationLog(String). Send message", "notification_log_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getNotificationLogWithCount")
    @Consumes("application/x-www-form-urlencoded")
    public Response getNotificationLogWithCount(@FormParam("json") String json) {     
        answer = CreateMessage(RestNotificationLogEnum.Notification_log_select_with_count.ordinal(), json, "getNotificationLogWithCount(String). Send message", "notification_log_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("addNotificationLog")
    @Consumes("application/x-www-form-urlencoded")
    public Response addNotificationLog(@FormParam("json") String json) {
        answer = CreateMessage(RestNotificationLogEnum.Notification_log_add.ordinal(), json, "addNotificationLog(String). Send message", "notification_log_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteNotificationLog")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteNotificationLog(@FormParam("json") String json) {
        answer = CreateMessage(RestNotificationLogEnum.Notification_log_delete.ordinal(), json, "deleteNotificationLog(String). Send message", "notification_log_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateNotificationLog")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateNotificationLog(@FormParam("json") String json) {
        answer = CreateMessage(RestNotificationLogEnum.Notification_log_update.ordinal(), json, "updateNotificationLog(String). Send message", "notification_log_queue");
        return InitMessage(answer);
    }
    
    @Override
    public ResData CreateMessage(int code, String json, String ErrorMessage, String name_microservice) {
        try {        
            JsonParser parser = new JsonParser();
            JsonObject json_obj = parser.parse(json).getAsJsonObject();
            String user_id = checkFildString(json_obj, "user_id");
            
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            return message_manager.call(builder.create().toJson(new RestObject(json_obj, user_id, code)), name_microservice);
        }
        catch (Exception ex) {
            return log_manager.CallError(ex.toString(), ErrorMessage, name_class);
        }
    } 
}
