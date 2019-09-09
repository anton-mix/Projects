package RestObject;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonParser;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class RessiveClassConverter implements JsonSerializer<RessiveClass>, JsonDeserializer<RessiveClass>, InterfaceConverter  {
    @Override
    public JsonElement serialize(RessiveClass src, Type type, JsonSerializationContext jsc) {
        JsonObject session_object = new com.google.gson.JsonObject();
        session_object.addProperty("json_message", src.getJson_message().toString());
        session_object.addProperty("code", src.getCode());
        return session_object;
    }
    
    @Override
    public RessiveClass deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        String json_message1 = checkFildString(object, "json_message"); 
        int code = checkFildInt(object, "code"); 
        JsonParser parser = new JsonParser();
        JsonObject json_message = parser.parse(json_message1).getAsJsonObject();

        return new RessiveClass(json_message, code);
    }
}

