package services;

import Managers.LogManager;
import Managers.MessageManager;
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
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.core.Response;

@Path("ManagerContacts")
public class ManagerContactsRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerContactsRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerContactsRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerContactRest");
        
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            log_manager.CallError(ex.toString(), "ManagerContactsRest. Init", name_class);
        }
    }
    
    
    @POST
    @Path("getContact")
    @Consumes("application/x-www-form-urlencoded")
    public Response getContact(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Contact_select_one.ordinal(), json, "getContact(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getContacts")
    @Consumes("application/x-www-form-urlencoded")
    public Response getContacts(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Contact_select_view.ordinal(), json, "getContacts(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getCountContacts")
    @Consumes("application/x-www-form-urlencoded")
    public Response getCountContacts(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Contact_select_count.ordinal(), json, "getCountContacts(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getContactsWithCount")
    public Response getContactsWithCount(@FormParam("json") String json) {     
        answer = CreateMessage(RestMailingEnum.Contact_select_with_count.ordinal(), json, "getContactsWithCount(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("addContact")
    @Consumes("application/x-www-form-urlencoded")
    public Response addContact(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Contact_add.ordinal(), json, "addContact(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteContact")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteContact(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Contact_delete.ordinal(), json, "deleteContact(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateContact")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateContact(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Contact_update.ordinal(), json, "updateContact(String). Send message", "mailing_session_queue");
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