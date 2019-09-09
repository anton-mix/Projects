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
import java.util.ArrayList;

public class TypeGroupContactsConverter implements JsonSerializer<TypeGroupContacts>, JsonDeserializer<TypeGroupContacts>, InterfaceConverter {
        @Override
    public JsonElement serialize(TypeGroupContacts src, Type type, JsonSerializationContext jsc) {
        JsonObject group_contact = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            group_contact.addProperty(key[i], values[i]);
        
        JsonArray groups_contacts = new JsonArray();
        for(int i = 0; i < src.getGroups_contacts().size();i++) {
            JsonObject group_contacts = new JsonObject();
            group_contacts.addProperty("id", src.getGroups_contacts().get(i).getId());
            group_contacts.addProperty("group_name", src.getGroups_contacts().get(i).getGroup_name());
            groups_contacts.add(group_contacts);
        }
        
        group_contact.add("groups_contacts", groups_contacts);
        
        JsonArray sender = new JsonArray();
        for(int i = 0; i < src.getSender().size();i++) {
            JsonObject email = new JsonObject();
            email.addProperty("id", src.getSender().get(i).getId());
            email.addProperty("email", src.getSender().get(i).getEmail());
            sender.add(email);
        }
        
        group_contact.add("sender", sender);
        
        return group_contact;
    }

    @Override
    public TypeGroupContacts deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        String name = checkFildString(object, "name");
        int sender_id = checkFildInt(object, "sender_id");
        int default_group_id = checkFildInt(object, "default_group_id");
        String telegram_bot_token = checkFildString(object, "telegram_bot_token");
        String vk_group_token = checkFildString(object, "vk_group_token");
        String vk_group_name = checkFildString(object, "vk_group_name");
        int vk_group_id = checkFildInt(object, "vk_group_id");
        
        ArrayList<ShortGroupContacts> groups_contacts = new ArrayList<ShortGroupContacts>();
        JsonArray arr = checkFildArrayJson(object,"groups_contacts");
        for(int i = 0; i < arr.size(); i++) {
            groups_contacts.add(new ShortGroupContacts(((JsonObject)arr.get(i)).get("id").getAsInt(), ((JsonObject)arr.get(i)).get("group_name").getAsString()));
        }
        
        return new TypeGroupContacts(id, user_id, name, sender_id, groups_contacts, default_group_id, telegram_bot_token, vk_group_token, vk_group_name, vk_group_id);
    }
}