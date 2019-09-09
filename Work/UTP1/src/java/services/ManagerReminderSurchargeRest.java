package services;

import Managers.LogManager;
import Managers.MessageManager;
import MyEnum.RestShopEnum;
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

@Path("ManagerReminderSurcharges")
public class ManagerReminderSurchargeRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerReminderSurchargeRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerReminderSurchargeRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerReminderSurchargeRest");
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            answer = log_manager.CallError(ex.toString(), "ManagerReminderSurchargeRest(). Init", name_class);
        } 
    }
    
    @POST
    @Path("getReminderSurcharges")
    @Consumes("application/x-www-form-urlencoded")
    public Response getReminderSurcharge(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_surcharge_select.ordinal(), json, "getReminderSurcharge(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getCountReminderSurcharges")
    @Consumes("application/x-www-form-urlencoded")
    public Response getCountReminderSurcharge(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_surcharge_select_count.ordinal(), json, "getCountReminderSurcharge(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getReminderSurchargesWithCount")
    @Consumes("application/x-www-form-urlencoded")
    public Response getReminderSurchargeWithCount(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_surcharge_select_with_count.ordinal(), json, "getReminderSurchargeWithCount(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("addReminderSurcharge")
    @Consumes("application/x-www-form-urlencoded")
    public Response addReminderSurcharge(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_surcharge_add.ordinal(), json, "addReminderSurcharge(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteReminderSurcharge")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteReminderSurcharge(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_surcharge_delete.ordinal(), json, "deleteReminderSurcharge(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateReminderSurcharge")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateReminderSurcharge(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_surcharge_update.ordinal(), json, "updateReminderSurcharge(String). Send message", "shop_session_queue");
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
