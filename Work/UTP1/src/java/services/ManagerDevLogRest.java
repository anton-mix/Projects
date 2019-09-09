package services;

import Managers.LogManager;
import Managers.MessageManager;
import MyEnum.DevLogEnum;
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

@Path("ManagerCallBack")
public class ManagerDevLogRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerDevLogRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerDevLogRest() {
        global_variables = new GlobalVariables();
          
        log_manager = new LogManager(global_variables, "ManagerDevLogRest");
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            log_manager.CallError(ex.toString(), "ManagerDevLogRest. Init", name_class);
        }
    }
    
    @POST
    @Path("addCallback")
    @Consumes("application/x-www-form-urlencoded")
    public Response addCallback(@FormParam("json") String json) {
        answer = CreateMessage(DevLogEnum.Write_log.ordinal(), json, "addCallback(String). Send message", "dev_log_queue");
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
