package QueueClass;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class UserAuthModelConverter implements JsonSerializer<UserAuthModel>, JsonDeserializer<UserAuthModel>, InterfaceConverter{
    
    @Override
    public JsonElement serialize(UserAuthModel src, Type type, JsonSerializationContext jsc) {
        JsonObject user = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            user.addProperty(key[i], values[i]);
      
        return user;
    }
    
    @Override
    public UserAuthModel deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        String name = checkFildString(object, "name");
        String username = checkFildString(object, "username");
        String password = checkFildString(object, "password");
        String session_id = checkFildString(object, "session_id");
        String email = checkFildString(object, "email");
        String phone = checkFildString(object, "phone");
        
        return new UserAuthModel(name, username, password, session_id, email, phone);
    }
}
