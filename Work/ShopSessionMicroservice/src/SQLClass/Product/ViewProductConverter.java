package SQLClass.Product;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class ViewProductConverter implements JsonSerializer<ViewProduct>, JsonDeserializer<ViewProduct>, InterfaceConverter  {

    @Override
    public JsonElement serialize(ViewProduct src, Type type, JsonSerializationContext jsc) {
        JsonObject viewProduct = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            viewProduct.addProperty(key[i], values[i]);
        
        return viewProduct;
    }
    
    @Override
    public ViewProduct deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object,"id");
        int user_id = checkFildInt(object,"user_id");
        String product_type = checkFildString(object,"product_type");
        int payment_type = checkFildInt(object,"payment_type");
        String product_name = checkFildString(object,"product_name");
        String identifier = checkFildString(object,"identifier");
        double value = checkFildDouble(object,"value");
        return new ViewProduct(id, user_id, product_type, payment_type, product_name, identifier, value);
    }
}
