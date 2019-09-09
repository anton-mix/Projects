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
import java.sql.Timestamp;

public class ProductModelConverter implements JsonSerializer<ProductModel>, JsonDeserializer<ProductModel>, InterfaceConverter {

    @Override
    public JsonElement serialize(ProductModel src, Type type, JsonSerializationContext jsc) {
        JsonObject product = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            product.addProperty(key[i], values[i]);
        return product;
    }

    @Override
    public ProductModel deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object,"id");
        int user_id = checkFildInt(object,"user_id");
        int group_id = checkFildInt(object,"group_id");
        String product_type = checkFildString(object,"product_type");
        int payment_type = checkFildInt(object,"payment_type");
        String order_page = checkFildString(object,"order_page");
        String product_name = checkFildString(object,"product_name");
        String identifier = checkFildString(object,"identifier");
        String image = checkFildString(object,"image");
        String product_description = checkFildString(object,"product_description");
        String groups_after_prepayment = checkFildString(object,"groups_after_prepayment");
        String message_topic_after_prepayment = checkFildString(object,"message_topic_after_prepayment");
        String message_type_after_prepayment = checkFildString(object,"message_type_after_prepayment");
        String message_after_prepayment = checkFildString(object,"message_after_prepayment");
        String groups_after_payment = checkFildString(object,"groups_after_payment");
        String thanks_page = checkFildString(object,"thanks_page");
        String message_topic_after_payment = checkFildString(object,"message_topic_after_payment");
        String message_type_after_payment = checkFildString(object,"message_type_after_payment");
        String message_after_payment = checkFildString(object,"message_after_payment");
        String getting_page = checkFildString(object,"getting_page");
        String pincodes = checkFildString(object,"pincodes");
        String vendor_code = checkFildString(object,"vendor_code");
        String company = checkFildString(object,"company");
        String country = checkFildString(object,"country");
        String ordering = checkFildString(object,"ordering");
        String bill_type = checkFildString(object,"bill_type");
        String employee_instruction = checkFildString(object,"employee_instruction");
        String personal_template = checkFildString(object,"personal_template");
        String some_payments = checkFildString(object,"some_payments");
        String show_in_catalogue = checkFildString(object,"show_in_catalogue");
        String tax = checkFildString(object,"tax");
        String show_for_partners = checkFildString(object,"show_for_partners");
        String personal_topic = checkFildString(object,"personal_topic");
        String personal_text = checkFildString(object,"personal_text");
        double value = checkFildDouble(object,"value");
        double prepayment_min = checkFildDouble(object,"prepayment_min");
        double amount_of_expenses = checkFildDouble(object,"amount_of_expenses");
        double percents_of_expenses = checkFildDouble(object,"percents_of_expenses");
        double comission = checkFildDouble(object,"comission");
        double comission_percents = checkFildDouble(object,"comission_percents");
        double employee_reward = checkFildDouble(object,"employee_reward");
        double employee_reward_percents = checkFildDouble(object,"employee_reward_percents");
        Timestamp bill_date = checkFildTimestamp(object,"bill_date");
        int reminder_payment_id = checkFildInt(object, "reminder_payment_id");
        int reminder_surcharge_id = checkFildInt(object, "reminder_surcharge_id");

        return new ProductModel (id, user_id, group_id, product_type, payment_type, order_page, product_name, identifier,value, some_payments, prepayment_min,
                image, show_in_catalogue, product_description, amount_of_expenses, percents_of_expenses, groups_after_prepayment, message_topic_after_prepayment,
                message_type_after_prepayment, message_after_prepayment, groups_after_payment, thanks_page, message_topic_after_payment,
                message_type_after_payment, message_after_payment, getting_page, pincodes, tax, vendor_code, company, country, ordering,
                bill_type, bill_date, comission, comission_percents, show_for_partners, employee_reward, employee_reward_percents, employee_instruction,
                personal_template, personal_topic, personal_text, reminder_payment_id, reminder_surcharge_id);

    }
}
