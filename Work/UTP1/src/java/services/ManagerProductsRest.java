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

@Path("ManagerProducts")
public class ManagerProductsRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerProductsRest";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;
    Gson gson;
    
    public ManagerProductsRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerProductsRest");
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            answer = log_manager.CallError(ex.toString(), "ManagerProductsRest(). Init", name_class);
        } 
    }
    
    @GET
    @Path("getProducts")
    public Response getProducts() {
        answer = CreateMessage(RestShopEnum.Product_select_all.ordinal(), "", "getProducts(). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getProducts")
    @Consumes("application/x-www-form-urlencoded")
    public Response getProducts(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Product_select_view.ordinal(), json, "getProducts(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getProduct")
    @Consumes("application/x-www-form-urlencoded")
    public Response getProduct(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Product_select_one.ordinal(), json, "getProduct(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getShortProducts")
    @Consumes("application/x-www-form-urlencoded")
    public Response getShortProducts(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Product_select_short.ordinal(), json, "getShortProducts(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getCountProducts")
    @Consumes("application/x-www-form-urlencoded")
    public Response getCountProducts(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Product_select_count.ordinal(), json, "getCountProducts(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("getProductsWithCount")
    @Consumes("application/x-www-form-urlencoded")
    public Response getProductsWithCount(@FormParam("json") String json) {     
        answer = CreateMessage(RestShopEnum.Product_select_with_count.ordinal(), json, "getCountProducts(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("addProduct")
    @Consumes("application/x-www-form-urlencoded")
    public Response addProduct(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Product_add.ordinal(), json, "addProduct(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("deleteProduct")
    @Consumes("application/x-www-form-urlencoded")
    public Response deleteProduct(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Product_delete.ordinal(), json, "deleteProduct(String). Send message", "shop_session_queue");
        return InitMessage(answer);
    }
    
    @POST
    @Path("updateProduct")
    @Consumes("application/x-www-form-urlencoded")
    public Response updateProduct(@FormParam("json") String json) {
        answer = CreateMessage(RestShopEnum.Product_update.ordinal(), json, "updateProduct(String). Send message", "shop_session_queue");
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