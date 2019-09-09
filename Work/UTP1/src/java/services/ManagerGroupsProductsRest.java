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

@Path("ManagerGroupsProducts")
public class ManagerGroupsProductsRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerGroupsProductsRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerGroupsProductsRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerGroupsProductsRest");
        try {
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            answer = log_manager.CallError(ex.toString(), "ManagerGroupsProductsRest. Init", name_class);
        } 
    }

    @POST
    @Path("getGroupsProducts")
    @Consumes("application/x-www-form-urlencoded")
    public Response getGroupsProducts(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Group_product_select.ordinal(), json, "getGroupsProducts(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getCountGroupsProducts")
    @Consumes("application/x-www-form-urlencoded")
    public Response getCountGroupsProducts(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Group_product_select_count.ordinal(), json, "getCountGroupsProducts(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("addGroupProducts")
    @Consumes("application/x-www-form-urlencoded")
    public Response addGroupProducts(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Group_product_add.ordinal(), json, "addGroupProducts(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteGroupProducts")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteGroupProducts(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Group_product_delete.ordinal(), json, "deleteGroupProducts(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateGroupProducts")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateGroupProducts(@FormParam("json") String json) { 
        answer = CreateMessage(RestShopEnum.Group_product_update.ordinal(), json, "updateGroupProducts(String). Send message", "shop_session_queue");
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