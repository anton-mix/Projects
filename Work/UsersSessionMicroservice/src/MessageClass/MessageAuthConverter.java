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

public class MessageAuthConverter implements JsonSerializer<MessageAuth>, JsonDeserializer<MessageAuth>, InterfaceConverter{
    
    @Override
    public JsonElement serialize(MessageAuth src, Type type, JsonSerializationContext jsc) {
        JsonObject session_object = new JsonObject();
        session_object.addProperty("name", src.getName());
        session_object.addProperty("session_id", src.getSession_id());
        session_object.addProperty("toEmail", src.getToEmail());
        session_object.addProperty("fromEmail", src.getFromEmail());
        return session_object;
    }
    
    @Override
    public MessageAuth deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        String name = checkFildString(object, "name"); 
        String session_id = checkFildString(object, "session_id"); 
        String toEmail = checkFildString(object, "toEmail"); 
        String fromEmail = checkFildString(object, "fromEmail"); 
        
        return new MessageAuth(name, session_id, toEmail, fromEmail);
    }
}
