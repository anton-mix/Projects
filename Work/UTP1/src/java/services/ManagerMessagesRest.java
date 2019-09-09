package services;

import Managers.LogManager;
import Managers.MessageManager;
import MyEnum.MailMessageEnum;
import MyEnum.RestMailingEnum;
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
import javax.ws.rs.GET;
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.core.Response;

@Path("ManagerMessages")
public class ManagerMessagesRest implements RestInterface{
    ResData answer;
    String name_class = "ManagerMessagesRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerMessagesRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerMessagesRest");
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            log_manager.CallError(ex.toString(), "ManagerMessagesRest(). Init", name_class);
        }
    }

    @POST
    @Path("getMessage")
    @Consumes("application/x-www-form-urlencoded")
    public Response getMessage(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Messages_select_one.ordinal(), json, "getMessage(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getMessages")
    @Consumes("application/x-www-form-urlencoded")
    public Response getMessages(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Messages_select_view.ordinal(), json, "getMessages(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getCountMessage")
    @Consumes("application/x-www-form-urlencoded")
    public Response getCountMessage(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Messages_select_count.ordinal(), json, "getCountMessage(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getMessagesWithCount")
    public Response getMessagesWithCount(@FormParam("json") String json) {   
        answer = CreateMessage(RestMailingEnum.Messages_select_with_count.ordinal(), json, "getMessagesWithCount(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
   
    @POST
    @Path("addMessage")
    @Consumes("application/x-www-form-urlencoded")
    public Response addMessage(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Message_add.ordinal(), json, "addMessage(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteMessage")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteMessage(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Message_delete.ordinal(), json, "deleteMessage(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateMessage")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateMessage(@FormParam("json") String json) {  
        answer = CreateMessage(RestMailingEnum.Message_update.ordinal(), json, "updateMessage(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getHttpSettings")
    @Consumes("application/x-www-form-urlencoded")
    public Response getHttpSettings(@FormParam("json") String json) {  
        answer = CreateMessage(RestMailingEnum.Http_settings_select.ordinal(), json, "getHttpSettings(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("addHttpSettings")
    @Consumes("application/x-www-form-urlencoded")
    public Response addHttpSettings(@FormParam("json") String json) {  
        answer = CreateMessage(RestMailingEnum.Http_settings_add.ordinal(), json, "addHttpMessage(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteHttpSettings")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteHttpSettings(@FormParam("json") String json) {  
        answer = CreateMessage(RestMailingEnum.Http_settings_delete.ordinal(), json, "deleteHttpSettings(String). Send messages", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateHttpSettings")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateHttpSettings(@FormParam("json") String json) {  
        answer = CreateMessage(RestMailingEnum.Http_settings_update.ordinal(), json, "updateHttpSettings(String). Send messages", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("sendMessage")
    @Consumes("application/x-www-form-urlencoded")
    public Response sendMessage(@FormParam("json") String json) {  
        answer = CreateMessage(MailMessageEnum.SendOtherMessage.ordinal(), json, "sendMessage(String). Send messages",  "mail_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("sendTgMessage")
    @Consumes("application/x-www-form-urlencoded")
    public Response sendTgMessage(@FormParam("json") String json) {  
        answer = CreateMessage(1, json, "sendMessage(String). Send messages", "tg_queue");
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
