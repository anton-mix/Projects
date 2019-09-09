package MyInterface;

import SQLClass.ResData;
import com.google.gson.JsonObject;
import javax.ws.rs.core.Response;

public interface RestInterface {
   ResData CreateMessage(int code, String json, String ErrorMessage, String name_microservice);
  
   default String checkFildString(JsonObject object, String fild) {
        try {
            return object.get(fild).getAsString();
        }
        catch(NullPointerException e){
            return "";
        }
    }
   
   default Response InitMessage(ResData answer) {
       return Response.ok()
                .header("Access-Control-Allow-Origin", "*")
                .header("Access-Control-Allow-Methods", "GET, POST, DELETE, PUT")
                .entity(answer.getRes()).build();
   } 
}
