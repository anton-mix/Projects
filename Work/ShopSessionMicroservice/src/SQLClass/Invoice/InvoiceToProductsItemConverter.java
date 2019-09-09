package SQLClass.Invoice;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class InvoiceToProductsItemConverter implements JsonSerializer<InvoiceToProductsItem>, JsonDeserializer<InvoiceToProductsItem>, InterfaceConverter {

    @Override
    public JsonElement serialize(InvoiceToProductsItem src, Type type, JsonSerializationContext jsc) {
        JsonObject invoicesToProductsItem = new JsonObject();
        invoicesToProductsItem.addProperty("id", src.getId());
        invoicesToProductsItem.addProperty("user_id", src.getUser_id());
        invoicesToProductsItem.addProperty("invoice_id", src.getInvoice_id());
        invoicesToProductsItem.addProperty("product_id", src.getProduct_id());
        return invoicesToProductsItem;
    }
    
    @Override
    public InvoiceToProductsItem deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        int invoice_id = checkFildInt(object, "invoice_id");
        int product_id = checkFildInt(object, "product_id");
        return new InvoiceToProductsItem(id, user_id, invoice_id, product_id);
    }
}
