package SQLClass.GroupObject;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;


public class GroupObjectConverter implements JsonSerializer<GroupObject>, JsonDeserializer<GroupObject>, InterfaceConverter  {
    @Override
    public JsonElement serialize(GroupObject src, Type type, JsonSerializationContext jsc) {
        JsonObject session_object = new JsonObject();
        session_object.addProperty("id", src.getId());
        session_object.addProperty("name", src.getName());
        session_object.addProperty("screen_name", src.getScreen_name());
        session_object.addProperty("is_closed", src.getIs_closed());
        session_object.addProperty("type", src.getType());
        session_object.addProperty("is_admin", src.getIs_admin());
        session_object.addProperty("admin_level", src.getAdmin_level());
        session_object.addProperty("is_member", src.getIs_member());
        session_object.addProperty("is_advertiser", src.getIs_advertiser());
        session_object.addProperty("photo_50", src.getPhoto_50());
        session_object.addProperty("photo_100", src.getPhoto_100());
        session_object.addProperty("photo_200", src.getPhoto_200());
        return session_object;
    }
    
    @Override
    public GroupObject deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = object.get("id").getAsInt(); 
        String name = object.get("name").getAsString(); 
        String screen_name = object.get("screen_name").getAsString(); 
        int is_closed = object.get("is_closed").getAsInt(); 
        String typeg = object.get("type").getAsString(); 
        int is_admin = object.get("is_admin").getAsInt(); 
        int admin_level = object.get("admin_level").getAsInt(); 
        int is_member = object.get("is_member").getAsInt(); 
        int is_advertiser = object.get("is_advertiser").getAsInt(); 
        String photo_50 = object.get("photo_50").getAsString(); 
        String photo_100 = object.get("photo_100").getAsString(); 
        String photo_200 = object.get("photo_200").getAsString(); 
        return new GroupObject(id, name, screen_name, is_closed, typeg, is_admin,
                admin_level, is_member, is_advertiser, photo_50, photo_100, photo_200);
    }
}
