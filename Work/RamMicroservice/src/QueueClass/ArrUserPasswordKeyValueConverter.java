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

public class ArrUserPasswordKeyValueConverter implements JsonSerializer<ArrUserPasswordKeyValue>, JsonDeserializer<ArrUserPasswordKeyValue>, InterfaceConverter {
    @Override
    public JsonElement serialize(ArrUserPasswordKeyValue src, Type type, JsonSerializationContext jsc) {
        JsonObject arrUserPasswordKey = new JsonObject();
        arrUserPasswordKey.addProperty("date", src.getDate().toString());
        arrUserPasswordKey.addProperty("login", src.getLogin());
        return arrUserPasswordKey;
    }
    
    @Override
    public ArrUserPasswordKeyValue deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        Timestamp date = checkFildTimestamp(object, "date");
        String login = checkFildString(object, "login");
        return new ArrUserPasswordKeyValue(date, login);
    }
}
