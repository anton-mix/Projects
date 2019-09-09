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

public class ContactEmailConverter implements JsonSerializer<ContactEmail>, JsonDeserializer<ContactEmail>, InterfaceConverter{

    @Override
    public JsonElement serialize(ContactEmail src, Type type, JsonSerializationContext jsc) {
        JsonObject contact_email = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            contact_email.addProperty(key[i], values[i]);
        return contact_email;
    }

    @Override
    public ContactEmail deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id"); 
        int user_id = checkFildInt(object, "user_id"); 
        int contact_id = checkFildInt(object, "contact_id"); 
        String email = checkFildString(object, "email"); 
        return new ContactEmail(id, user_id, contact_id, email);
    }
    
}
