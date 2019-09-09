package SendAndResivClass;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class ResultConverter implements JsonSerializer<Result>, JsonDeserializer<Result>, InterfaceConverter{

    @Override
    public JsonElement serialize(Result src, Type type, JsonSerializationContext jsc) {
        JsonObject result = new JsonObject();
        result.addProperty("result", src.getResult());
        return result;
    }

    @Override
    public Result deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        String res = checkFildString(object, "result");
        return new Result(res);
    }
}
