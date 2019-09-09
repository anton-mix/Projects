package SQLClass.Social;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class VkGroupsConverter  implements JsonSerializer<VkGroups>, JsonDeserializer<VkGroups>, InterfaceConverter {
    
    @Override
    public JsonElement serialize(VkGroups src, Type type, JsonSerializationContext jsc) {
        JsonObject userSystemInfo = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            userSystemInfo.addProperty(key[i], values[i]);
        
        return userSystemInfo;
    }
    
    @Override
    public VkGroups deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        int type_id = checkFildInt(object, "type_id");
        String group_name = checkFildString(object, "group_name");
        String group_token = checkFildString(object, "group_token");
        int group_id = checkFildInt(object, "group_id");
        int enabled = checkFildInt(object, "enabled");
        int banned = checkFildInt(object, "banned");

        return new VkGroups(id, user_id, type_id, group_name, group_token, group_id, enabled, banned);
    }
}
