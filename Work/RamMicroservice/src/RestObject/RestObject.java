package RestObject;

import MyEnum.RestRamEnum;
import com.google.gson.JsonObject;

public class RestObject {
    JsonObject json_message;
    String session_key = "";
    RestRamEnum code;
    
    public RestObject(JsonObject json_message, String session_key, RestRamEnum code){
        this.json_message = json_message;
        this.session_key = session_key;
        this.code = code;
    }

    public JsonObject getJson_message() {
        return json_message;
    }

    public void setJson_message(JsonObject json_message) {
        this.json_message = json_message;
    }

    public String getSession_key() {
        return session_key;
    }

    public void setSession_key(String session_key) {
        this.session_key = session_key;
    }

    public RestRamEnum getCode() {
        return code;
    }

    public void setCode(RestRamEnum code) {
        this.code = code;
    }    
    
    @Override
    public String toString(){
        return "Json message: " + json_message.toString() + " Session key: " + session_key + " Code: " + code.name();
    }
}

