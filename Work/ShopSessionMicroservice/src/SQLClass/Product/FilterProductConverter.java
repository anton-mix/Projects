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
import java.util.List;

public class FilterProductConverter implements JsonSerializer<FilterProduct>, JsonDeserializer<FilterProduct>, InterfaceConverter {

    @Override
    public JsonElement serialize(FilterProduct src, Type type, JsonSerializationContext jsc) {
        JsonObject filterProduct = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            filterProduct.addProperty(key[i], values[i]);
        
        JsonArray categories = new JsonArray();
        for(int i = 0; i < src.getCategories_id().size();i++) {
            JsonObject categorie = new JsonObject();
            categorie.addProperty("id", src.getCategories_id().get(i).getId());
            categorie.addProperty("group_name",src.getCategories_id().get(i).getGroup_name());
            categories.add(categorie);
        }
        filterProduct.add("categories_id", categories);
        return filterProduct;
    }

    @Override
    public FilterProduct deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {  
        JsonObject object = json.getAsJsonObject();
        int user_id = checkFildInt(object,"user_id");
        String show_in_catalogue = checkFildString(object,"show_in_catalogue");
        String query = checkFildString(object,"query");
        int limit = checkFildInt(object,"limit");
        int offset = checkFildInt(object,"offset");
        List<ShortGroupProducts> categories_id = new ArrayList<ShortGroupProducts>();
        JsonArray arr = checkFildArrayJson(object,"categories_id");
        for(int i = 0; i < arr.size();i++) {
            categories_id.add(new ShortGroupProducts(((JsonObject)arr.get(i)).get("id").getAsInt(), ((JsonObject)arr.get(i)).get("group_name").getAsString()));
        }
        return new FilterProduct(user_id, show_in_catalogue, query, categories_id, limit, offset);
    }
    
}
