package SQLClass.VKAuthClass;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;


public class VKAuthClassConverter implements JsonSerializer<VKAuthClass>, JsonDeserializer<VKAuthClass>, InterfaceConverter{

    @Override
    public JsonElement serialize(VKAuthClass src, Type type, JsonSerializationContext jsc) {
        JsonObject vk_auth = new JsonObject();
        vk_auth.addProperty("token", src.getToken());
        vk_auth.addProperty("id", src.getId());
        vk_auth.addProperty("user_id", src.getUser_id());
        return vk_auth;
    }

    @Override
    public VKAuthClass deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        String token = checkFildString(object, "token");
        int user_id = checkFildInt(object, "user_id");
        return new VKAuthClass(id, user_id, token);
    }
    
}
