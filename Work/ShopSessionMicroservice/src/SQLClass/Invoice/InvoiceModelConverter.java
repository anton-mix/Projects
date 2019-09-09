package SQLClass.Invoice;

import MyInterface.InterfaceConverter;
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
import java.util.ArrayList;

public class InvoiceModelConverter implements JsonSerializer<InvoiceModel>, JsonDeserializer<InvoiceModel>, InterfaceConverter {

    @Override
    public JsonElement serialize(InvoiceModel src, Type type, JsonSerializationContext jsc) {
        JsonObject invoice = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            invoice.addProperty(key[i], values[i]);
        
        JsonArray products = new JsonArray();
        for(int i = 0; i < src.products.size();i++) {
            JsonObject product = new JsonObject();
            product.addProperty("id", src.products.get(i).getId());
            product.addProperty("value", src.products.get(i).getValue());
            product.addProperty("product_name",src.products.get(i).getProduct_name());
            products.add(product);
        }
        invoice.add("products", products);
        
        return invoice;
    }
    
    @Override
    public InvoiceModel deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object,"id");
        String email = checkFildString(object,"email");
        String name = checkFildString(object,"name");
        String phone = checkFildString(object,"phone");
        String comment = checkFildString(object,"comment");
        String telegram = checkFildString(object,"telegram");
        String vk = checkFildString(object,"vk");
        String facebook = checkFildString(object,"facebook");
        String postcode = checkFildString(object,"postcode");
        String country = checkFildString(object,"country");
        String region = checkFildString(object,"region");
        String city = checkFildString(object,"city");
        String address = checkFildString(object,"address");
        String channel = checkFildString(object,"channel");
        String source = checkFildString(object,"source");
        String campaign = checkFildString(object,"campaign");
        String advertising_ads = checkFildString(object,"advertising_ads");
        String advertising_keys = checkFildString(object,"advertising_keys");
        int user_id = checkFildInt(object,"user_id");
        String customer_time = checkFildString(object,"customer_time");
        String account_validity_period_bool = checkFildString(object,"account_validity_period_bool");
        double amount = checkFildDouble(object,"amount");
        Timestamp account_validity_period = checkFildTimestamp(object,"account_validity_period");
        Timestamp date_create = checkFildTimestamp(object,"date_create");
        String payment_status = checkFildString(object,"payment_status");
        String ip_address = checkFildString(object, "ip_address");
        ArrayList<ShortProduct> products = new ArrayList<ShortProduct>();
        JsonArray arr = checkFildArrayJson(object,"products");
        for(int i = 0; i < arr.size();i++) {
            products.add(new ShortProduct(
                    ((JsonObject)arr.get(i)).get("id").getAsInt(), 
                    ((JsonObject)arr.get(i)).get("product_name").getAsString()));
        }

        return new InvoiceModel(id, user_id, email, name, phone, amount, customer_time, account_validity_period,
            account_validity_period_bool, date_create, comment, telegram, vk, facebook, postcode,
            country, region, city, address, channel, source, campaign, advertising_ads, 
            advertising_keys, payment_status, products, ip_address);
    }
}
