
package SQLClass.VkError;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class VkErrorConverter implements JsonSerializer<VkError>, JsonDeserializer<VkError>, InterfaceConverter  {
    @Override
    public JsonElement serialize(VkError src, Type type, JsonSerializationContext jsc) {
        JsonObject session_object = new JsonObject();
        session_object.addProperty("error", src.getErr());
        return session_object;
    }
    
    @Override
    public VkError deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        String json_message = checkFildJson(object, "error").toString();
        return new VkError(json_message);
    }
}
