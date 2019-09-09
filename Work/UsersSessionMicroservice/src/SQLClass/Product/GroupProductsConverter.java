package SQLClass.Product;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonArray;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;
import java.util.ArrayList;

public class GroupProductsConverter implements JsonSerializer<GroupProducts>, JsonDeserializer<GroupProducts>, InterfaceConverter {

    @Override
    public JsonElement serialize(GroupProducts src, Type type, JsonSerializationContext jsc) {
        JsonObject groupProducts = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            groupProducts.addProperty(key[i], values[i]);
        
        JsonArray products = new JsonArray();
        for(int i = 0; i < src.products.size();i++) {
            JsonObject product = new JsonObject();
            product.addProperty("id", src.products.get(i).getId());
            product.addProperty("value",src.products.get(i).getValue());
            product.addProperty("product_name",src.products.get(i).getProduct_name());
            products.add(product);
        }
        groupProducts.add("products", products);
        return groupProducts;
    }

    @Override
    public GroupProducts deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object,"id");
        int user_id = checkFildInt(object,"user_id");
        String group_name = checkFildString(object,"group_name");
        ArrayList<ShortProduct> products = new ArrayList<ShortProduct>();
        JsonArray arr = checkFildArrayJson(object,"products");
        for(int i = 0; i < arr.size();i++) {
            products.add(new ShortProduct(((JsonObject)arr.get(i)).get("id").getAsInt(),((JsonObject)arr.get(i)).get("product_name").getAsString()));
        }
        return new GroupProducts(id, user_id, group_name, products);
    }
    
}
