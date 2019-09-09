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

public class GroupContactsConverter  implements JsonSerializer<GroupContacts>, JsonDeserializer<GroupContacts>, InterfaceConverter {
    
    @Override
    public JsonElement serialize(GroupContacts src, Type type, JsonSerializationContext jsc) {
        JsonObject group_contact = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            group_contact.addProperty(key[i], values[i]);
        
        JsonArray emails = new JsonArray();
        for(int i = 0; i < src.sender.size();i++) {
            JsonObject email = new JsonObject();
            email.addProperty("id", src.sender.get(i).getId());
            email.addProperty("email", src.sender.get(i).getEmail());
            emails.add(email);
        }
        group_contact.add("sender", emails);
        
        return group_contact;
    }

    @Override
    public GroupContacts deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        String type_name = checkFildString(object, "type_name");
        int type_id = checkFildInt(object, "type_id");
        String group_name = checkFildString(object, "group_name"); 
        int amount_contacts = checkFildInt(object, "amount_contacts");
        int activated_contacts = checkFildInt(object, "activated_contacts");
        String activation_text = checkFildString(object, "activation_text");
        String unsubscribe_url = checkFildString(object, "unsubscribe_url");
        int sender_id = checkFildInt(object, "sender_id");
        
        return new GroupContacts(id, user_id, group_name, type_id, type_name, amount_contacts, activated_contacts, activation_text, unsubscribe_url, sender_id);
    }
}
