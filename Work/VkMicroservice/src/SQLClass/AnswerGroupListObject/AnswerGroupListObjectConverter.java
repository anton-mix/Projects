package SQLClass.AnswerGroupListObject;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonArray;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;


public class AnswerGroupListObjectConverter implements JsonSerializer<AnswerGroupListObject>, JsonDeserializer<AnswerGroupListObject>, InterfaceConverter  {
    @Override
    public JsonElement serialize(AnswerGroupListObject src, Type type, JsonSerializationContext jsc) {
        JsonObject session_object = new JsonObject();
        session_object.addProperty("json_message", src.getJson_message().toString());
        session_object.addProperty("count", src.getCount());
        return session_object;
    }
    
    @Override
    public AnswerGroupListObject deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int count = object.get("count").getAsInt(); 
        JsonArray json_message = object.get("items").getAsJsonArray();
        return new AnswerGroupListObject(count, json_message);
    }
}
