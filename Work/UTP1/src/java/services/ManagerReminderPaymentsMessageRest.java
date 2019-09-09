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

@Path("ManagerReminderPaymentMails")
public class ManagerReminderPaymentsMessageRest implements RestInterface{
    ResData answer;
    String name_class = "ManagerReminderPaymentsMessageRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerReminderPaymentsMessageRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerReminderPaymentsMessageRest");
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            answer = log_manager.CallError(ex.toString(), "ManagerReminderPaymentsMessageRest(). Init", name_class);
        } 
    }
    
    @POST
    @Path("getReminderPaymentMails")
    @Consumes("application/x-www-form-urlencoded")
    public Response getReminderPaymentMails(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_payments_messages_select.ordinal(), json, "getReminderPaymentMails(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("addReminderPaymentMail")
    @Consumes("application/x-www-form-urlencoded")
    public Response addReminderPaymentsMail(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_payments_messages_add.ordinal(), json, "addReminderPaymentsMail(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteReminderPaymentMail")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteReminderPaymentsMail(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_payments_messages_delete.ordinal(), json, "deleteReminderPaymentsMail(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateReminderPaymentMail")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateReminderPaymentsMail(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_payments_messages_update.ordinal(), json, "updateReminderPaymentsMail(String). Send message", "shop_session_queue");
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
