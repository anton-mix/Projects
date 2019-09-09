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

@Path("ManagerReminderSurchargeMails")
public class ManagerReminderSurchargeMessageRest implements RestInterface {

    ResData answer;
    String name_class = "ManagerReminderSurchargeMessageRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;

    public ManagerReminderSurchargeMessageRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerReminderSurchargeMessageRest");

        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();

            message_manager = new MessageManager(global_variables);
        } catch (Exception ex) {
            answer = log_manager.CallError(ex.toString(), "ManagerReminderSurchargeMessageRest(). Init", name_class);
        }
    }

    @POST
    @Path("getReminderSurchargeMail")
    @Consumes("application/x-www-form-urlencoded")
    public Response getReminderSurchargeMail(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_surcharge_messages_select.ordinal(), json, "getReminderSurchargeMail(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }

    @POST
    @Path("addReminderSurchargeMail")
    @Consumes("application/x-www-form-urlencoded")
    public Response addReminderSurchargeMail(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_surcharge_messages_add.ordinal(), json, "addReminderSurchargeMail(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }

    @POST
    @Path("deleteReminderSurchargeMail")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteReminderSurchargeMail(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_surcharge_messages_delete.ordinal(), json, "deleteReminderSurchargeMail(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }

    @POST
    @Path("updateReminderSurchargeMail")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateReminderSurchargeMail(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Reminder_surcharge_messages_update.ordinal(), json, "updateReminderSurchargeMail(String). Send message", "shop_session_queue");
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
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), ErrorMessage, name_class);
        }
    }
}
