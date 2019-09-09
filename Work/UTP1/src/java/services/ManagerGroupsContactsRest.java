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

@Path("ManagerGroupsContacts")
public class ManagerGroupsContactsRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerGroupsContactsRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerGroupsContactsRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerGroupsContactsRest");
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            answer = log_manager.CallError(ex.toString(), "ManagerGroupsContactsRest. Init", name_class);
        } 
    }
    
    public Gson initGson(){
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
        gson = builder.create();
        return gson;
    }
    
    @POST
    @Path("getGroupsContacts")
    @Consumes("application/x-www-form-urlencoded")
    public Response getGroupsContacts(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Group_contact_select_view.ordinal(), json, "getGroupsContacts(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getGroupContacts")
    @Consumes("application/x-www-form-urlencoded")
    public Response getGroupContacts(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Group_contact_select_one.ordinal(), json, "getGroupContacts(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("addGroupContacts")
    @Consumes("application/x-www-form-urlencoded")
    public Response addGroupContacts(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Group_contact_add.ordinal(), json, "addGroupContacts(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteGroupContacts")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteGroupContacts(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Group_contact_delete.ordinal(), json, "deleteGroupContacts(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateGroupContacts")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateGroupContacts(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Group_contact_update.ordinal(), json, "updateGroupContacts(String). Send message", "mailing_session_queue");
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
