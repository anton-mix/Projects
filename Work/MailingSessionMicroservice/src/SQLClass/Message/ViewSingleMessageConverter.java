
package SQLClass.Message;

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

public class ViewSingleMessageConverter implements JsonSerializer<ViewSingleMessage>, JsonDeserializer<ViewSingleMessage>, InterfaceConverter{
    @Override
    public JsonElement serialize(ViewSingleMessage src, Type type, JsonSerializationContext jsc) {
        JsonObject view_message = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            view_message.addProperty(key[i], values[i]);
        return view_message;
    }
    
    @Override
    public ViewSingleMessage deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id"); 
        int user_id = checkFildInt(object, "user_id");
        String message_head = checkFildString(object, "message_head");
        Timestamp date_send = checkFildTimestamp(object, "date_send");
        String view_date = date_send.toString().split(" ", 2)[0];
        String view_time = date_send.toString().split(" ", 2)[1];    
        String amount_receivers = checkFildString(object, "amount_receivers");
        String amount_opening = checkFildString(object, "amount_opening");
        String amount_visiters = checkFildString(object, "amount_visiters");
        String unknown_field = checkFildString(object, "unknown_field");
        
        return new ViewSingleMessage(id, user_id, message_head, date_send, view_date, view_time, 
            amount_receivers, amount_opening, amount_visiters, unknown_field);
    }
    
    
}
