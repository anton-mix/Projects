package SQLClass.ReminderMaill;

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

public class ReminderMailConverter implements JsonSerializer<ReminderMail>, JsonDeserializer<ReminderMail>, InterfaceConverter{
    @Override
    public JsonElement serialize(ReminderMail src, Type type, JsonSerializationContext jsc) {
        JsonObject reminderMessageView = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            reminderMessageView.addProperty(key[i], values[i]);
        return reminderMessageView;
    }

    @Override
    public ReminderMail deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object,"id");
        int user_id = checkFildInt(object,"user_id");
        int contact_id = checkFildInt(object,"contact_id");
        int maill_tamplate_id = checkFildInt(object,"maill_tamplate_id");
        int product_id = checkFildInt(object,"product_id");
        Timestamp date_send = checkFildTimestamp(object,"date_send");
        return new ReminderMail(id, user_id, contact_id, maill_tamplate_id, product_id, date_send);
    }
}
