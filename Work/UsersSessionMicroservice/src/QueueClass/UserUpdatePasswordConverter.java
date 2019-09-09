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

public class UserUpdatePasswordConverter implements JsonSerializer<UserUpdatePassword>, JsonDeserializer<UserUpdatePassword>, InterfaceConverter {
    @Override
    public JsonElement serialize(UserUpdatePassword src, Type type, JsonSerializationContext jsc) {
        JsonObject userUpdatePassword = new JsonObject();
        userUpdatePassword.addProperty("user_id", src.getUser_id());
        userUpdatePassword.addProperty("old_password",src.getOld_password());
        userUpdatePassword.addProperty("new_password",src.getNew_password());
        return userUpdatePassword;
    }
    
    @Override
    public UserUpdatePassword deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int user_id = checkFildInt(object, "user_id");
        String old_password = checkFildString(object, "old_password");
        String new_password = checkFildString(object, "new_password");     
        return new UserUpdatePassword(user_id, old_password, new_password);
    }
}
