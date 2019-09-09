/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SQLClass.NotificationLog;

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

public class NotificationLogModelConverter implements JsonSerializer<NotificationLogModel>, JsonDeserializer<NotificationLogModel>, InterfaceConverter {
    
    @Override
    public JsonElement serialize(NotificationLogModel src, Type type, JsonSerializationContext jsc) {
        JsonObject contact_email = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            contact_email.addProperty(key[i], values[i]);
        return contact_email;
    }

    @Override
    public NotificationLogModel deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id"); 
        int user_id = checkFildInt(object, "user_id"); 
        Timestamp date = checkFildTimestamp(object, "date"); 
        String content = checkFildString(object, "content"); 
        return new NotificationLogModel(id, user_id, date, content);
    }
}
