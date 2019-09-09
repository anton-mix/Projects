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

public class SessionArrValueConverter implements JsonSerializer<SessionArrValue>, JsonDeserializer<SessionArrValue>, InterfaceConverter {
    @Override
    public JsonElement serialize(SessionArrValue src, Type type, JsonSerializationContext jsc) {
        JsonObject sessionArrValue = new JsonObject();
        sessionArrValue.addProperty("id", src.getId());
        sessionArrValue.addProperty("date", src.getDate().toString());
        return sessionArrValue;
    }
    
    @Override
    public SessionArrValue deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        Timestamp date = checkFildTimestamp(object, "date");
        return new SessionArrValue(id, date);
    }
}
