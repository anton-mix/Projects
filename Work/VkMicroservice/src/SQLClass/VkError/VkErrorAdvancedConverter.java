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

public class VkErrorAdvancedConverter implements JsonSerializer<VkErrorAdvanced>, JsonDeserializer<VkErrorAdvanced>, InterfaceConverter  {
    @Override
    public JsonElement serialize(VkErrorAdvanced src, Type type, JsonSerializationContext jsc) {
        JsonObject session_object = new JsonObject();
        session_object.addProperty("error_code", src.getCode());
        session_object.addProperty("error_msg", src.getError_msg());
        return session_object;
    }
    
    @Override
    public VkErrorAdvanced deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int code = object.get("error_code").getAsInt(); 
        String error_msg = object.get("error_msg").getAsString(); 
        return new VkErrorAdvanced(code, error_msg);
    }
}