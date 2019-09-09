package SQLClass.Contact;

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

public class ViewGroupContactsConverter implements JsonSerializer<ViewGroupContacts>, JsonDeserializer<ViewGroupContacts>, InterfaceConverter {

    @Override
    public JsonElement serialize(ViewGroupContacts src, Type type, JsonSerializationContext jsc) {
        JsonObject group_contact = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            group_contact.addProperty(key[i], values[i]);
        return group_contact;
    }

    @Override
    public ViewGroupContacts deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        int type_id = checkFildInt(object, "type_id");
        String type_name = checkFildString(object, "type_name");
        String group_name = checkFildString(object, "group_name"); 
        int amount_contacts = checkFildInt(object, "amount_contacts");
        int activated_contacts = checkFildInt(object, "activated_contacts");
        
        return new ViewGroupContacts(id, user_id, type_id, group_name, type_name, amount_contacts, activated_contacts);
    }
    
}
