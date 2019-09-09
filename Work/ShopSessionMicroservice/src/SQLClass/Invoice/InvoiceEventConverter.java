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
import java.sql.Timestamp;

public class InvoiceEventConverter implements JsonSerializer<InvoiceEvent>, JsonDeserializer<InvoiceEvent>, InterfaceConverter {

    @Override
    public JsonElement serialize(InvoiceEvent src, Type type, JsonSerializationContext jsc) {
        JsonObject invoiceEvent = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            invoiceEvent.addProperty(key[i], values[i]);
        return invoiceEvent;
    }

    @Override
    public InvoiceEvent deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id"); 
        int user_id = checkFildInt(object, "user_id"); 
        int invoice_id = checkFildInt(object, "invoice_id"); 
        int typef = checkFildInt(object, "type"); 
        String is_system = checkFildString(object,"is_system");
        String title = checkFildString(object, "title");
        Timestamp date = checkFildTimestamp(object, "date");
        String owner = checkFildString(object, "owner"); 
        String is_call_back = checkFildString(object, "is_call_back");
        String call_back = checkFildString(object, "call_back");
        String comment = checkFildString(object, "comment");
        String view_date = date.toString().split(" ", 2)[0];
        String view_time = date.toString().split(" ", 2)[1];
        
        return new InvoiceEvent(id, user_id, invoice_id, typef, is_system, title, date, owner, is_call_back, call_back, comment, view_date, view_time);
    }   
}
