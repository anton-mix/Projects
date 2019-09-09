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

@Path("ManagerHttpSettings")
public class ManagerHttpSettingsRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerHttpSettingsRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerHttpSettingsRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerHttpSettingsRest");
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            answer = log_manager.CallError(ex.toString(), "ManagerHttpSettingsRest(). Init", name_class);
        } 
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
    public Response AddHttpSettings(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Http_settings_add.ordinal(), json, "addHttpSettings(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteHttpSettings")
    @Consumes("application/x-www-form-urlencoded")
    public Response DeleteHttpSettings(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Http_settings_delete.ordinal(), json, "deleteHttpSettings(String). Send message", "mailing_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateHttpSettings")
    @Consumes("application/x-www-form-urlencoded")
    public Response UpdateHttpSettings(@FormParam("json") String json) {
        answer = CreateMessage(RestMailingEnum.Http_settings_update.ordinal(), json, "updateHttpSettings(String). Send message", "mailing_session_queue");
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

