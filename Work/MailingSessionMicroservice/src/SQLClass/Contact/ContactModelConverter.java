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
import java.sql.Timestamp;
import java.util.ArrayList;

public class ContactModelConverter implements JsonSerializer<ContactModel>, JsonDeserializer<ContactModel>, InterfaceConverter {

    @Override
    public JsonElement serialize(ContactModel src, Type type, JsonSerializationContext jsc) {
        JsonObject contact = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            contact.addProperty(key[i], values[i]);
        
        JsonArray groupsContacts = new JsonArray();
        for(int i = 0; i < src.groups_contacts.size();i++) {
            JsonObject groupContacts = new JsonObject();
            groupContacts.addProperty("id", src.groups_contacts.get(i).getId());
            groupContacts.addProperty("group_name",src.groups_contacts.get(i).getGroup_name());
            groupsContacts.add(groupContacts);
        }
        contact.add("groups_contacts", groupsContacts);
        
        JsonArray socials = new JsonArray();
        for(int i = 0; i < src.social.size();i++) {
            JsonObject social = new JsonObject();
            social.addProperty("value", src.social.get(i));
            socials.add(social);
        }
        contact.add("social", socials);
        
        JsonArray emails = new JsonArray();
        for(int i = 0; i < src.emails.size();i++) {
            JsonObject email = new JsonObject();
            email.addProperty("value", src.emails.get(i));
            emails.add(email);
        }
        contact.add("emails", emails);
        
        JsonArray phones = new JsonArray();
        for(int i = 0; i < src.phones.size();i++) {
            JsonObject phone = new JsonObject();
            phone.addProperty("value", src.phones.get(i));
            phones.add(phone);
        }
        contact.add("phones", phones);
       
        return contact;
    }

    @Override
    public ContactModel deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id"); 
        int user_id = checkFildInt(object, "user_id"); 
        String surname = checkFildString(object, "surname"); 
        String name = checkFildString(object, "name"); 
        String patronymic = checkFildString(object, "patronymic"); 
        String city = checkFildString(object, "city"); 
        String client_info = checkFildString(object, "client_info"); 
        String tags = checkFildString(object, "tags"); 
        String contact_status = checkFildString(object, "contact_status"); 
        Timestamp date_create = checkFildTimestamp(object,"date_create");
        ArrayList<ShortGroupContacts> groups_contacts = new ArrayList<ShortGroupContacts>();
        JsonArray arr = checkFildArrayJson(object,"groups_contacts");
        for(int i = 0; i < arr.size();i++) {
            groups_contacts.add(new ShortGroupContacts(((JsonObject)arr.get(i)).get("id").getAsInt(), ((JsonObject)arr.get(i)).get("group_name").getAsString()));
        }
        
        ArrayList<String> social = new ArrayList<String>();
        JsonArray arr2 = checkFildArrayJson(object,"social");
        for(int i = 0; i < arr2.size();i++) {
            social.add(((JsonObject)arr2.get(i)).get("value").getAsString());
        }
        
        ArrayList<String> emails = new ArrayList<String>();
        JsonArray arr3 = checkFildArrayJson(object,"emails");
        for(int i = 0; i < arr3.size();i++) {
            emails.add(((JsonObject)arr3.get(i)).get("value").getAsString());
        }
        
        ArrayList<String> phones = new ArrayList<String>();
        JsonArray arr4 = checkFildArrayJson(object,"phones");
        for(int i = 0; i < arr4.size();i++) {
            phones.add(((JsonObject)arr4.get(i)).get("value").getAsString());
        }
        
        return new ContactModel(id, user_id, surname, name, patronymic, city, client_info, tags, contact_status, date_create, groups_contacts,
            social, emails, phones);
    }
    
}
