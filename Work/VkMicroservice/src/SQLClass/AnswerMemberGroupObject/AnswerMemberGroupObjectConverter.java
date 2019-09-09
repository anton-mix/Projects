package SQLClass.AnswerMemberGroupObject;

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
import java.util.ArrayList;


public class AnswerMemberGroupObjectConverter implements JsonSerializer<AnswerMemberGroupObject>, JsonDeserializer<AnswerMemberGroupObject>, InterfaceConverter  {
    @Override
    public JsonElement serialize(AnswerMemberGroupObject src, Type type, JsonSerializationContext jsc) {
        JsonObject session_object = new JsonObject();
        session_object.addProperty("count", src.getCount());
        session_object.addProperty("ids", src.getIds().toString());
        return session_object;
    }
    
    @Override
    public AnswerMemberGroupObject deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int count = object.get("count").getAsInt(); 
        JsonArray jsonArray = json.getAsJsonObject().getAsJsonArray("items");
        ArrayList<String> ids = new ArrayList<String>();
        for(int i = 0; i < jsonArray.size();i ++) {
            ids.add(jsonArray.get(i).getAsString());
        }
        return new AnswerMemberGroupObject(count, ids);
    }
}