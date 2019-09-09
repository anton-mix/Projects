package SQLClass;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class SQLResultConverter implements JsonSerializer<SQLResult>, JsonDeserializer<SQLResult>, InterfaceConverter{

    @Override
    public JsonElement serialize(SQLResult src, Type type, JsonSerializationContext jsc) {
        JsonObject sQLResult = new JsonObject();
        sQLResult.addProperty("result", src.getResult());
        sQLResult.addProperty("id", src.getId());
        sQLResult.addProperty("user_id", src.getUser_id());
        return sQLResult;
    }

    @Override
    public SQLResult deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        String res = checkFildString(object, "result");
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        return new SQLResult(res,id,user_id);
    }
    
}
