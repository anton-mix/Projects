package SQLClass.Contact;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class ContactPhoneConverter implements JsonSerializer<ContactPhone>, JsonDeserializer<ContactPhone>, InterfaceConverter{
    @Override
    public JsonElement serialize(ContactPhone src, Type type, JsonSerializationContext jsc) {
        JsonObject contact_email = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            contact_email.addProperty(key[i], values[i]);
        return contact_email;
    }

    @Override
    public ContactPhone deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id"); 
        int user_id = checkFildInt(object, "user_id"); 
        int contact_id = checkFildInt(object, "contact_id"); 
        String phone = checkFildString(object, "phone"); 
        return new ContactPhone(id, user_id, contact_id, phone);
    }
}
