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

public class VkContactsConverter implements JsonSerializer<VkContacts>, JsonDeserializer<VkContacts>, InterfaceConverter {
    
    @Override
    public JsonElement serialize(VkContacts src, Type type, JsonSerializationContext jsc) {
        JsonObject userSystemInfo = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            userSystemInfo.addProperty(key[i], values[i]);
        
        return userSystemInfo;
    }
    
    @Override
    public VkContacts deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        int contact_id  = checkFildInt(object, "contact_id");
        String first_name = checkFildString(object, "first_name");
        String last_name = checkFildString(object, "last_name");
        int is_subscribed = checkFildInt(object, "is_subscribed");
        
        return new VkContacts(id, user_id, contact_id, first_name, last_name, is_subscribed);
    }
}
