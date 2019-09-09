package SQLClass.VKSendingClass;

import MyEnum.RestVkEnum;
import MyInterface.InterfaceConverter;
import RestObject.RestObject;
import com.google.gson.JsonArray;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonParser;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;
import java.util.ArrayList;

public class VKSendingClassConverter implements JsonSerializer<VKSendingClass>, JsonDeserializer<VKSendingClass>, InterfaceConverter  {
    @Override
    public JsonElement serialize(VKSendingClass src, Type type, JsonSerializationContext jsc) {
        JsonObject session_object = new JsonObject();
        String message_text = src.getMessage_text();
        String u_tok = src.getUser_token();
        String g_tok = src.getGroup_token();
        int owner_id = src.getOwner_id();
        String group_id = src.getGroup_id();
        session_object.addProperty("message_text", message_text);
        session_object.addProperty("user_token", u_tok);
        session_object.addProperty("group_token", g_tok);
        session_object.addProperty("owner_id", owner_id);
        session_object.addProperty("group_id", group_id);
        JsonArray user_ids = new JsonArray();
        for(int i = 0; i < src.getUser_ids().size(); i++) {
            JsonObject user_id = new JsonObject();
            user_id.addProperty("user_id", src.getUser_ids().get(i));
            user_ids.add(user_id);
        }
        session_object.add("user_ids", user_ids);
        return session_object;
    }
    
    @Override
    public VKSendingClass deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        String message_text = checkFildString(object, "message_text"); 
        String u_tok = checkFildString(object,"user_token");
        String g_tok = checkFildString(object,"group_token");
        int owner_id = checkFildInt(object,"owner_id");
        String group_id = checkFildString(object,"group_id");
        JsonArray arr2 = checkFildArrayJson(object,"user_ids");
        
        ArrayList<Integer> user_ids = new ArrayList<Integer>();
        for(int i = 0; i < arr2.size();i++) {
            user_ids.add(((JsonObject)arr2.get(i)).get("user_id").getAsInt());
        }
        return new VKSendingClass(u_tok, g_tok, message_text, user_ids, owner_id, group_id);
    }
}
