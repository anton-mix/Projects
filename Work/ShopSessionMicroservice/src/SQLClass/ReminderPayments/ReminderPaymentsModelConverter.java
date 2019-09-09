package SQLClass.ReminderPayments;

import MyInterface.InterfaceConverter;
import SQLClass.Product.ShortProduct;
import SQLClass.ReminderMaill.ReminderMessageView;
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

public class ReminderPaymentsModelConverter implements JsonSerializer<ReminderPaymentsModel>, JsonDeserializer<ReminderPaymentsModel>, InterfaceConverter {
    @Override
    public JsonElement serialize(ReminderPaymentsModel src, Type type, JsonSerializationContext jsc) {
        JsonObject reminderPayments = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            reminderPayments.addProperty(key[i], values[i]);
        
        JsonArray products = new JsonArray();
        for(int i = 0; i < src.getProducts().size();i++) {
            JsonObject product = new JsonObject();
            product.addProperty("id", src.getProducts().get(i).getId());
            product.addProperty("product_name", src.getProducts().get(i).getProduct_name());
            products.add(product);
        }
        reminderPayments.add("products", products);
        
        JsonArray messages = new JsonArray();
        for(int i = 0; i < src.getMessages().size();i++) {
            JsonObject message = new JsonObject();
            message.addProperty("id", src.getMessages().get(i).getId());
            message.addProperty("message_topic", src.getMessages().get(i).getMessage_topic());
            message.addProperty("date_send", src.getMessages().get(i).getDate_send().toString());
            messages.add(message);
        }
        reminderPayments.add("messages", messages);
        
        return reminderPayments;
    }

    @Override
    public ReminderPaymentsModel deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object,"id");
        int user_id = checkFildInt(object,"user_id");
        String name = checkFildString(object,"name");
        int count = checkFildInt(object,"count");
        ArrayList<ShortProduct> products = new ArrayList<ShortProduct>();
        JsonArray arr = checkFildArrayJson(object,"products");
        for(int i = 0; i < arr.size();i++) {
            products.add(new ShortProduct(checkFildInt(arr.get(i).getAsJsonObject(),"id"),  checkFildString(arr.get(i).getAsJsonObject(),"product_name")));
        }
        
        ArrayList<ReminderMessageView> messages = new ArrayList<ReminderMessageView>();
        JsonArray arr2 = checkFildArrayJson(object,"messages");
        for(int i = 0; i < arr2.size();i++) {
            messages.add(new ReminderMessageView(checkFildInt(arr2.get(i).getAsJsonObject(),"id"),  checkFildString(arr2.get(i).getAsJsonObject(),"message_topic"), checkFildTimestamp(arr2.get(i).getAsJsonObject(), "date_send")));
        }
        return new ReminderPaymentsModel (id, user_id, name, count, products, messages);
    }
}
