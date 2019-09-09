package services;

import Managers.LogManager;
import Managers.MessageManager;
import MyEnum.RestMailingEnum;
import MyEnum.RestVkEnum;
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

@Path("ManagerTypesGroupsContacts")
public class ManagerTypesGroupsContactsRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerTypesGroupsContactsRests";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerTypesGroupsContactsRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerTypesGroupsContactsRest");
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            answer = log_manager.CallError(ex.toString(), "ManagerTypesGroupsContactsRest(). Init", name_class);
        } 
    }
    
    @POST
    @Path("getTypesGroupsContacts")
    @Consumes("application/x-www-form-urlencoded")
    public Response getvk(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Types_groups_contacts_select.ordinal(), json, "getTypesGroupsContacts(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("addTypeGroupContacts")
    @Consumes("application/x-www-form-urlencoded")
    public Response addTypeGroupContacts(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Types_groups_contacts_add.ordinal(), json, "addTypeGroupContacts(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteTypeGroupContacts")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteTypeGroupContacts(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Types_groups_contacts_delete.ordinal(), json, "deleteTypeGroupContacts(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateTypeGroupContacts")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateTypeGroupContacts(@FormParam("json") String json) { 
        answer = CreateMessage(RestMailingEnum.Types_groups_contacts_update.ordinal(), json, "updateTypeGroupContacts(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getVkGroups")
    @Consumes("application/x-www-form-urlencoded")
    public Response getVkGroup(@FormParam("json") String json) { 
        answer = CreateMessage(RestVkEnum.get_user_groups.ordinal(), json, "getVkGroup(String). Send message", "vk_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("SendVkMessage")
    @Consumes("application/x-www-form-urlencoded")
    public Response sendVkGroup(@FormParam("json") String json) { 
        answer = CreateMessage(RestVkEnum.send_messages_by_group.ordinal(), json, "SendVkMessage(String). Send message", "vk_queue");
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