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
import java.sql.Timestamp;

public class UserArrValueConverter implements JsonSerializer<UserArrValue>, JsonDeserializer<UserArrValue>, InterfaceConverter {
    @Override
    public JsonElement serialize(UserArrValue src, Type type, JsonSerializationContext jsc) {
        JsonObject userArrValue = new JsonObject();
        userArrValue.addProperty("date", src.getDate().toString());
        JsonObject user = new JsonObject();
        user.addProperty("name", src.getUser().getName());
        user.addProperty("username", src.getUser().getUsername());
        user.addProperty("password", src.getUser().getPassword());
        user.addProperty("session_id", src.getUser().getSession_id());
        user.addProperty("email", src.getUser().getEmail());
        user.addProperty("phone", src.getUser().getPhone());
        userArrValue.add("user", user);
        return userArrValue;
    }
    
    @Override
    public UserArrValue deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        Timestamp date = checkFildTimestamp(object, "date");
        JsonObject user = checkFildJson(object, "user");
        UserAuthModel buffuser = new UserAuthModel(checkFildString(user, "name"), checkFildString(user, "username") ,checkFildString(user, "password"),
        checkFildString(user, "session_id"), checkFildString(user, "email"), checkFildString(user, "phone"));
        return new UserArrValue(date, buffuser);
    }
}
