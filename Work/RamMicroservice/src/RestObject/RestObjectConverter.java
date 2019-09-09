package RestObject;

import MyEnum.RestRamEnum;
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

public class RestObjectConverter implements JsonSerializer<RestObject>, JsonDeserializer<RestObject>, InterfaceConverter  {
    @Override
    public JsonElement serialize(RestObject src, Type type, JsonSerializationContext jsc) {
        JsonObject session_object = new JsonObject();
        session_object.addProperty("json_message", src.getJson_message().toString());
        session_object.addProperty("session_key", src.getSession_key());
        session_object.addProperty("code", src.getCode().ordinal());
        return session_object;
    }
    
    @Override
    public RestObject deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        String json_message1 = object.get("json_message").getAsString(); 
        JsonParser parser = new JsonParser();
        JsonObject json_message = parser.parse(json_message1).getAsJsonObject();
        String session_key = checkFildString(object, "session_key");
        RestRamEnum code = RestRamEnum.values()[checkFildInt(object, "code")];
        return new RestObject(json_message, session_key, code);
    }
}
