package SQLClass.ReminderMaill;

import MyInterface.InterfaceConverter;
import SQLClass.Product.GroupProducts;
import SQLClass.Product.ShortProduct;
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

public class ReminderMessageViewConverter implements JsonSerializer<ReminderMessageView>, JsonDeserializer<ReminderMessageView>, InterfaceConverter {
    
    @Override
    public JsonElement serialize(ReminderMessageView src, Type type, JsonSerializationContext jsc) {
        JsonObject reminderMessageView = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            reminderMessageView.addProperty(key[i], values[i]);
        return reminderMessageView;
    }

    @Override
    public ReminderMessageView deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object,"id");
        String message_topic = checkFildString(object,"message_topic");
        Timestamp date_send = checkFildTimestamp(object,"date_send");
        return new ReminderMessageView(id, message_topic, date_send);
    }
}
