package SQLClass.Social;

import MyInterface.InterfaceConverter;
import SQLClass.User.UserSystemInfo;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class TelegramBotsConverter implements JsonSerializer<TelegramBots>, JsonDeserializer<TelegramBots>, InterfaceConverter {
    @Override
    public JsonElement serialize(TelegramBots src, Type type, JsonSerializationContext jsc) {
        JsonObject userSystemInfo = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            userSystemInfo.addProperty(key[i], values[i]);
        
        return userSystemInfo;
    }
    
    @Override
    public TelegramBots deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        int type_id = checkFildInt(object, "type_id");
        String bot_token = checkFildString(object, "bot_token");
        String bot_name = checkFildString(object, "bot_name");
        int enabled = checkFildInt(object, "enabled");
        int banned = checkFildInt(object, "banned");
        return new TelegramBots(id, user_id, type_id, bot_token, bot_name, enabled, banned);
    }
}
