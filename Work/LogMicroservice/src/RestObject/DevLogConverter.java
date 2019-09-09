package RestObject;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonParser;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;
import java.sql.Timestamp;

public class DevLogConverter  implements JsonSerializer<DevLog>, JsonDeserializer<DevLog>, InterfaceConverter {
    @Override
    public JsonElement serialize(DevLog src, Type type, JsonSerializationContext jsc) {
        JsonObject dev_log_object = new JsonObject(); 
        dev_log_object.addProperty("id", src.getId());
        dev_log_object.addProperty("user_id", src.getUser_id());
        dev_log_object.addProperty("date", src.getDate().toString());
        dev_log_object.addProperty("theme", src.getTheme());
        dev_log_object.addProperty("content", src.getContent());
        return dev_log_object;
    }
    
    @Override
    public DevLog deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        Timestamp date = checkFildTimestamp(object, "date");;
        String theme = checkFildString(object, "theme");
        String content = checkFildString(object, "content");
        return new DevLog(id, user_id, date, theme, content);
    }
}
