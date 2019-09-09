package MyInterface;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import java.sql.Timestamp;

public interface InterfaceConverter {
    
    default String checkFildString(JsonObject object, String fild) {
        try {
            return object.get(fild).getAsString();
        }
        catch(NullPointerException e){
            return "";
        }
    }
    
    default int checkFildInt(JsonObject object, String fild) {
        try {            
            return object.get(fild).getAsInt();
        }
        catch(NullPointerException e){
            return -1;
        }
    }
    
    default Double checkFildDouble(JsonObject object, String fild) {
        try {
            return object.get(fild).getAsDouble();
        }
        catch(NullPointerException e){
            return -1.0;
        }
    }
    
    default JsonArray checkFildArrayJson(JsonObject object, String fild) {
        try {
            return object.get(fild).getAsJsonArray();
        }
        catch(NullPointerException e){
            return new JsonArray();
        }
    }
    
    default Timestamp checkFildTimestamp(JsonObject object, String fild) {
        try {
            return Timestamp.valueOf(object.get(fild).getAsString());
        }
        catch(Exception e) {
            return Timestamp.valueOf("1970-01-02 0:0:0");
        }
    }
    
    default JsonObject checkFildJson(JsonObject object, String fild) {
        try {
            return object.get(fild).getAsJsonObject();
        }
        catch(NullPointerException e){
            return new JsonObject();
        }
    }
}
