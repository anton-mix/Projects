package RestObject;

import com.google.gson.JsonObject;

public class RessiveClass {
    JsonObject json_message;
    int code;
    
    public RessiveClass(JsonObject json_message, int code){
        this.json_message = json_message;
        this.code = code;
    }

    public JsonObject getJson_message() {
        return json_message;
    }

    public void setJson_message(JsonObject json_message) {
        this.json_message = json_message;
    }

    public int getCode() {
        return code;
    }

    public void setCode(int code) {
        this.code = code;
    }
    
    @Override
    public String toString() {
        return "Json message: " + json_message.toString() + " Code:" + code;
    }
}