package MessageClass;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class MessageChangePasswordConverter implements JsonSerializer<MessageChangePassword>, JsonDeserializer<MessageChangePassword>, InterfaceConverter {
        @Override
    public JsonElement serialize(MessageChangePassword src, Type type, JsonSerializationContext jsc) {
        JsonObject session_object = new JsonObject();      
        session_object.addProperty("passwordKey", src.getPasswordKey());
        session_object.addProperty("toEmail", src.getToEmail());
        session_object.addProperty("fromEmail", src.getFromEmail());
        return session_object;
    }
    
    @Override
    public MessageChangePassword deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        String passwordKey = checkFildString(object, "passwordKey"); 
        String toEmail = checkFildString(object, "toEmail"); 
        String fromEmail = checkFildString(object, "fromEmail"); 
        return new MessageChangePassword(passwordKey, toEmail, fromEmail);
    }
}
