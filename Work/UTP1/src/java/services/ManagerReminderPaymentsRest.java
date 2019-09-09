package services;

import Managers.LogManager;
import Managers.MessageManager;
import MyEnum.RestNotificationLogEnum;
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

@Path("ManagerReminderPayments")
public class ManagerReminderPaymentsRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerReminderPaymentsRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerReminderPaymentsRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerReminderPaymentsRest");
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            answer = log_manager.CallError(ex.toString(), "ManagerReminderPaymentsRest(). Init", name_class);
        } 
    }
    
    @POST
    @Path("getReminderPayments")
    @Consumes("application/x-www-form-urlencoded")
    public Response getReminderPayments(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_payments_select.ordinal(), json, "getReminderPayments(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getCountReminderPayments")
    @Consumes("application/x-www-form-urlencoded")
    public Response getCountReminderPayments(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_payments_select_count.ordinal(), json, "getCountReminderPayments(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getReminderPaymentsWithCount")
    @Consumes("application/x-www-form-urlencoded")
    public Response getReminderPaymentsWithCount(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_payments_select_with_count.ordinal(), json, "getReminderPaymentsWithCount(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("addReminderPayment")
    @Consumes("application/x-www-form-urlencoded")
    public Response addReminderPayments(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_payments_add.ordinal(), json, "addReminderPayments(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteReminderPayment")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteReminderPayments(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_payments_delete.ordinal(), json, "deleteReminderPayments(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateReminderPayment")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateReminderPayments(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_payments_update.ordinal(), json, "updateReminderPayments(String). Send message", "shop_session_queue");
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
