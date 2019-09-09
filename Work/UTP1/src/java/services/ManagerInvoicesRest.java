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
import javax.ws.rs.GET;
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.core.Response;

@Path("ManagerInvoices")
public class ManagerInvoicesRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerInvoicesRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerInvoicesRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerInvoicesRest");
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            log_manager.CallError(ex.toString(), "ManagerInvoicesRest(). Init", name_class);
        }
    }
    
    @GET
    @Path("getInvoice")
    public Response getInvoice() {     
        answer = CreateMessage(RestShopEnum.Invoice_select_all.ordinal(), "",  "getInvoice2(). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getInvoice")
    @Consumes("application/x-www-form-urlencoded")
    public Response getInvoice(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Invoice_select_one.ordinal(), json, "getInvoice(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getInvoices")
    @Consumes("application/x-www-form-urlencoded")
    public Response getInvoices(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Invoice_select_view.ordinal(), json, "getInvoices(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getCountInvoice")
    @Consumes("application/x-www-form-urlencoded")
    public Response getCountInvoice(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Invoice_select_count.ordinal(), json, "getCountInvoice(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getInvoicesWithCount")
    public Response getInvoicesWithCount(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Invoice_select_with_count.ordinal(), json, "getInvoicesWithCount(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("addInvoice")
    @Consumes("application/x-www-form-urlencoded")
    public Response addInvoice(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Invoice_add.ordinal(), json, "addInvoice(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteInvoice")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteInvoice(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Invoice_delete.ordinal(), json, "deleteInvoice(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateInvoice")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateInvoice(@FormParam("json") String json) {  
        answer = CreateMessage(RestShopEnum.Invoice_update.ordinal(), json, "updateInvoice(String). Send message", "shop_session_queue");
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
            return message_manager.call(builder.create().toJson(new RestObject(json_obj, user_id, code)), "shop_session_queue");
        }
        catch (Exception ex) {
            return log_manager.CallError(ex.toString(), ErrorMessage, name_class);
        }
    } 
}
