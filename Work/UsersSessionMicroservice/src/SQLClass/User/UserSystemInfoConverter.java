package SQLClass.User;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class UserSystemInfoConverter implements JsonSerializer<UserSystemInfo>, JsonDeserializer<UserSystemInfo>, InterfaceConverter {
    @Override
    public JsonElement serialize(UserSystemInfo src, Type type, JsonSerializationContext jsc) {
        JsonObject userSystemInfo = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            userSystemInfo.addProperty(key[i], values[i]);
        
        return userSystemInfo;
    }
    
    @Override
    public UserSystemInfo deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        int group_products_id = checkFildInt(object, "group_products_id");
        int type_group_contacts = checkFildInt(object, "type_group_contacts_id");
        int reminder_payment_id_1 = checkFildInt(object, "reminder_payment_id_1");
        int reminder_payment_id_2 = checkFildInt(object, "reminder_payment_id_2");
        int reminder_payment_id_3 = checkFildInt(object, "reminder_payment_id_3");
        int reminder_surcharge_id_1 = checkFildInt(object, "reminder_payment_id_1");
        int reminder_surcharge_id_2 = checkFildInt(object, "reminder_payment_id_1");
        int reminder_surcharge_id_3 = checkFildInt(object, "reminder_payment_id_1");
        
        return new UserSystemInfo(id, user_id, group_products_id, type_group_contacts, 
                reminder_payment_id_1, reminder_payment_id_2, reminder_payment_id_3,
                reminder_surcharge_id_1,reminder_surcharge_id_2, reminder_surcharge_id_3);
    }
}
