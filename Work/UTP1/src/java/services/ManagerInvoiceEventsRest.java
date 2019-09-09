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

@Path("ManagerInvoiceEvent")
public class ManagerInvoiceEventsRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerHttpSettingsRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerInvoiceEventsRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerInvoiceEventsRest");
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            answer = log_manager.CallError(ex.toString(), "ManagerInvoiceEventsRest. Init", name_class);
        } 
    }

   
    @POST
    @Path("getInvoiceEvents")
    @Consumes("application/x-www-form-urlencoded")
    public Response getInvoiceEvents(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Invoice_event_select.ordinal(), json, "getInvoiceEvents(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getCountInvoiceEvents")
    @Consumes("application/x-www-form-urlencoded")
    public Response getCountInvoiceEvents(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Invoice_event_select_count.ordinal(), json, "getCountInvoiceEvents(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("addInvoiceEvent")
    @Consumes("application/x-www-form-urlencoded")
    public Response addInvoiceEvent(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Invoice_event_add.ordinal(), json, "addInvoiceEvent(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteInvoiceEvent")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteInvoiceEvent(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Invoice_event_delete.ordinal(), json, "deleteInvoiceEvent(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateInvoiceEvent")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateInvoiceEvent(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Invoice_event_update.ordinal(), json, "updateInvoiceEvent(String). Send message", "shop_session_queue");
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
