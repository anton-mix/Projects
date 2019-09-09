package SQLClass.Invoice;

import MyInterface.InterfaceConverter;
import SQLClass.Invoice.ViewInvoice;
import SQLClass.Product.ShortProduct;
import com.google.gson.JsonArray;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;
import java.sql.Timestamp;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

public class ViewInvoiceConverter implements JsonSerializer<ViewInvoice>, JsonDeserializer<ViewInvoice>, InterfaceConverter {
    @Override
    public JsonElement serialize(ViewInvoice src, Type type, JsonSerializationContext context) {
        JsonObject viewInvoice = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            viewInvoice.addProperty(key[i], values[i]);
        
        JsonArray products = new JsonArray();
        for(int i = 0; i < src.products.size();i++) {
            JsonObject product = new JsonObject();
            product.addProperty("id", src.products.get(i).getId());
            product.addProperty("value", src.products.get(i).getValue());
            product.addProperty("product_name",src.products.get(i).getProduct_name());
            products.add(product);
        }
        viewInvoice.add("products", products);
        
        return viewInvoice;
    }
    
    @Override
    public ViewInvoice deserialize(JsonElement json, Type type, JsonDeserializationContext context) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        String email = checkFildString(object, "email");
        String name = checkFildString(object, "name");
        String phone = checkFildString(object, "phone");
        int user_id = checkFildInt(object, "user_id");
        double amount = checkFildDouble(object, "amount");
        Timestamp date_create = checkFildTimestamp(object, "date_create");
        String payment_status = checkFildString(object, "payment_status");
        ArrayList<ShortProduct> products = new ArrayList<ShortProduct>();
        JsonArray arr = checkFildArrayJson(object, "products");
        for(int i = 0; i < arr.size();i++) {
            products.add(new ShortProduct(((JsonObject)arr.get(i)).get("id").getAsInt(), ((JsonObject)arr.get(i)).get("product_name").getAsString()));
        }
    
        try {
            return new ViewInvoice(id, user_id, email, name, phone, amount, date_create, payment_status, products);
        } catch (ParseException ex) {
            Logger.getLogger(ViewInvoiceConverter.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;
    }
}
