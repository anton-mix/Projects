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

public class ViewTypeGroupContactsConverter implements JsonSerializer<ViewTypeGroupContacts>, JsonDeserializer<ViewTypeGroupContacts>, InterfaceConverter {

    @Override
    public JsonElement serialize(ViewTypeGroupContacts src, Type type, JsonSerializationContext jsc) {
        JsonObject type_group_contacts = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            type_group_contacts.addProperty(key[i], values[i]);
        return type_group_contacts;
    }

    @Override
    public ViewTypeGroupContacts deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        String name = checkFildString(object, "name");
        
        return new ViewTypeGroupContacts(id,user_id, name);
    }
    
}