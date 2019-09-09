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

public class SimpleMessageConverter implements JsonSerializer<SimpleMessage>, JsonDeserializer<SimpleMessage>, InterfaceConverter {
        @Override
    public JsonElement serialize(SimpleMessage src, Type type, JsonSerializationContext jsc) {
        JsonObject session_object = new JsonObject();      
        session_object.addProperty("text", src.getText());
        session_object.addProperty("toEmail", src.getToEmail());
        session_object.addProperty("fromEmail", src.getFromEmail());
        return session_object;
    }
    
    @Override
    public SimpleMessage deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        String text = checkFildString(object, "text"); 
        String toEmail = checkFildString(object, "toEmail"); 
        String fromEmail = checkFildString(object, "fromEmail"); 
        return new SimpleMessage(text, toEmail, fromEmail);
    }
}
