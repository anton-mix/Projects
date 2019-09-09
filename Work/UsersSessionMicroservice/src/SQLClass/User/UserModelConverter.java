package SQLClass.User;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class UserModelConverter implements JsonSerializer<UserModel>, JsonDeserializer<UserModel>, InterfaceConverter {

    @Override
    public JsonElement serialize(UserModel src, Type type, JsonSerializationContext jsc) {
        JsonObject user = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            user.addProperty(key[i], values[i]);
        
        return user;
    }
    
    @Override
    public UserModel deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        String username = checkFildString(object, "username");
        String password = checkFildString(object, "password");
        String surname = checkFildString(object, "surname");
        String name = checkFildString(object, "name");
        String gpatronymic = checkFildString(object, "gpatronymic");
        String avatar_url = checkFildString(object, "avatar_url");
        String email = checkFildString(object, "email");
        String phone = checkFildString(object, "phone");
        
        return new UserModel(id, username, password, surname, name, gpatronymic, avatar_url, email, phone);
    }
}
