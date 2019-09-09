package SQLClass.Message;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonArray;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;
import java.sql.Timestamp;
import java.util.ArrayList;

public class FilterMessageConverter implements JsonSerializer<FilterMessage>, JsonDeserializer<FilterMessage>, InterfaceConverter {
    @Override
    public JsonElement serialize(FilterMessage src, Type type, JsonSerializationContext jsc) {
        JsonObject filterMessage = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            filterMessage.addProperty(key[i], values[i]);
           
        return filterMessage;
    }

    @Override
    public FilterMessage deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int user_id = checkFildInt(object, "user_id"); 
        int message_num = checkFildInt(object, "message_num"); 
        Timestamp date_send_begin = checkFildTimestamp(object, "date_send_begin");
        Timestamp date_send_end = checkFildTimestamp(object, "date_send_end");
        String message_topic = checkFildString(object, "message_topic");
        int limit = checkFildInt(object, "limit");
        int offset = checkFildInt(object, "offset");
        
        return new FilterMessage(user_id, message_num, date_send_begin, date_send_end, message_topic, limit, offset);
    }
}
