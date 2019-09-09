
package SQLClass.Message;

import MyInterface.InterfaceConverter;
import SQLClass.Contact.ViewGroupContacts;
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


public class SingleMessagesModelConverter implements JsonSerializer<SingleMessagesModel>, JsonDeserializer<SingleMessagesModel>, InterfaceConverter{
    @Override
    public JsonElement serialize(SingleMessagesModel src, Type type, JsonSerializationContext jsc) {
        JsonObject message = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            message.addProperty(key[i], values[i]);
                
        JsonArray sending_groups = new JsonArray();
        for(int i = 0; i < src.sending_groups.size();i++) {
            JsonObject sending_group = new JsonObject();
            sending_group.addProperty("id", src.sending_groups.get(i).getId());
            sending_group.addProperty("user_id", src.sending_groups.get(i).getUser_id());
            sending_group.addProperty("type_id", src.sending_groups.get(i).getType_id());
            sending_group.addProperty("type_name", src.sending_groups.get(i).getType_name());
            sending_group.addProperty("group_name", src.sending_groups.get(i).getGroup_name());
            sending_group.addProperty("amount_contacts", src.sending_groups.get(i).getAmount_contacts());
            sending_group.addProperty("activated_contacts", src.sending_groups.get(i).getActivated_contacts());
            sending_groups.add(sending_group);
        }
        message.add("sending_groups", sending_groups);
        
        JsonArray sending_invoices = new JsonArray();
        for(int i = 0; i < src.sending_invoices.size();i++) {
            JsonObject sending_invoice = new JsonObject();
            sending_invoice.addProperty("id", src.sending_invoices.get(i).getId());
            sending_invoice.addProperty("user_id", src.sending_invoices.get(i).getUser_id());
            sending_invoice.addProperty("email", src.sending_invoices.get(i).getEmail());
            sending_invoice.addProperty("name", src.sending_invoices.get(i).getName());
            sending_invoice.addProperty("phone", src.sending_invoices.get(i).getPhone());
            sending_invoice.addProperty("amount", src.sending_invoices.get(i).getAmount());
            sending_invoice.addProperty("date_create", src.sending_invoices.get(i).getDate_create().toString());
            sending_invoice.addProperty("payment_status", src.sending_invoices.get(i).getPayment_status());
            JsonArray products = new JsonArray();
            for (int j = 0; j < src.sending_invoices.get(j).getProducts().size(); j++){
                JsonObject product = new JsonObject();
                product.addProperty("id", src.sending_invoices.get(i).getProducts().get(j).getId());
                product.addProperty("value", src.sending_invoices.get(i).getProducts().get(j).getValue());
                product.addProperty("product_name", src.sending_invoices.get(i).getProducts().get(j).getProduct_name());
                products.add(product);
            }
            sending_invoice.add("products", products);
            sending_invoices.add(sending_invoice);
        }
        message.add("sending_invoices", sending_invoices);
        
        JsonArray sending_groups_exception = new JsonArray();
        for(int i = 0; i < src.sending_groups_exception.size();i++) {
            JsonObject one_sending_groups_exception = new JsonObject();
            one_sending_groups_exception.addProperty("id", src.sending_groups_exception.get(i).getId());
            one_sending_groups_exception.addProperty("user_id", src.sending_groups_exception.get(i).getUser_id());
            one_sending_groups_exception.addProperty("type_id", src.sending_groups_exception.get(i).getType_id());
            one_sending_groups_exception.addProperty("type_name", src.sending_groups_exception.get(i).getType_name());
            one_sending_groups_exception.addProperty("group_name", src.sending_groups_exception.get(i).getGroup_name());
            one_sending_groups_exception.addProperty("amount_contacts", src.sending_groups_exception.get(i).getAmount_contacts());
            one_sending_groups_exception.addProperty("activated_contacts", src.sending_groups_exception.get(i).getActivated_contacts());
            sending_groups_exception.add(one_sending_groups_exception);
        }
        message.add("sending_groups_exception", sending_groups_exception);
        
        JsonArray sending_invoices_exception = new JsonArray();
        for(int i = 0; i < src.sending_invoices_exception.size();i++) {
            JsonObject one_sending_invoices_exception = new JsonObject();
            one_sending_invoices_exception.addProperty("id", src.sending_invoices_exception.get(i).getId());
            one_sending_invoices_exception.addProperty("user_id", src.sending_invoices_exception.get(i).getUser_id());
            one_sending_invoices_exception.addProperty("email", src.sending_invoices_exception.get(i).getEmail());
            one_sending_invoices_exception.addProperty("name", src.sending_invoices_exception.get(i).getName());
            one_sending_invoices_exception.addProperty("phone", src.sending_invoices_exception.get(i).getPhone());
            one_sending_invoices_exception.addProperty("amount", src.sending_invoices_exception.get(i).getAmount());
            one_sending_invoices_exception.addProperty("date_create", src.sending_invoices_exception.get(i).getDate_create().toString());
            one_sending_invoices_exception.addProperty("payment_status", src.sending_invoices_exception.get(i).getPayment_status());
            
            JsonArray products = new JsonArray();
            for (int j = 0; j < src.sending_invoices_exception.get(j).getProducts().size(); j++){
                JsonObject product = new JsonObject();
                product.addProperty("id", src.sending_invoices_exception.get(i).getProducts().get(j).getId());
                product.addProperty("value", src.sending_invoices_exception.get(i).getProducts().get(j).getValue());
                product.addProperty("product_name", src.sending_invoices_exception.get(i).getProducts().get(j).getProduct_name());
                products.add(product);
            }
            one_sending_invoices_exception.add("products", products);
            sending_invoices_exception.add(one_sending_invoices_exception);
        }
        message.add("sending_invoices_exception", sending_invoices_exception);
        
        JsonArray message_links = new JsonArray();
        for(int i = 0; i < src.message_links.size();i++) {
            JsonObject message_link = new JsonObject();
            message_link.addProperty("id", src.message_links.get(i).getId());
            message_link.addProperty("user_id", src.message_links.get(i).getUser_id());
            message_link.addProperty("message_id", src.message_links.get(i).getMessage_id());
            message_link.addProperty("link_address", src.message_links.get(i).getLink_address());
            message_link.addProperty("link_name", src.message_links.get(i).getLink_name());
            message_link.addProperty("link_text_link", src.message_links.get(i).getLink_text_link());
            message_links.add(message_link);
        }
        message.add("message_links", message_links);
       
        JsonArray image_links = new JsonArray();
        for(int i = 0; i < src.image_links.size();i++) {
            JsonObject image_link = new JsonObject();
            image_link.addProperty("id", src.image_links.get(i).getId());
            image_link.addProperty("user_id", src.image_links.get(i).getUser_id());
            image_link.addProperty("message_id", src.image_links.get(i).getMessage_id());
            image_link.addProperty("image_address", src.image_links.get(i).getImage_address());
            image_link.addProperty("image_name", src.image_links.get(i).getImage_name());
            image_link.addProperty("image_text_link", src.image_links.get(i).getImage_text_link());
            image_links.add(image_link);
        }
        message.add("image_links", image_links);
        
        JsonArray messes = new JsonArray();
        for(int i = 0; i < src.message.size();i++) {
            JsonObject mess = new JsonObject();
            mess.addProperty("id", src.message.get(i).getId());
            mess.addProperty("user_id", src.message.get(i).getUser_id());
            mess.addProperty("type_id", src.message.get(i).getType_id());
            mess.addProperty("type_name", src.message.get(i).getType_name());
            mess.addProperty("group_name", src.message.get(i).getGroup_name());
            mess.addProperty("amount_contacts", src.message.get(i).getAmount_contacts());
            mess.addProperty("activated_contacts", src.message.get(i).getActivated_contacts());
            messes.add(mess);
        }
        message.add("message", messes);
        
        return message;
    }
    
    

    @Override
    public SingleMessagesModel deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id"); 
        int user_id = checkFildInt(object, "user_id");
        String status = checkFildString(object, "status");
        Timestamp date_send = checkFildTimestamp(object, "date_send");
        String view_date = date_send.toString().split(" ", 2)[0];
        String view_time = date_send.toString().split(" ", 2)[1];
        String time_zone = checkFildString(object, "time_zone");
        String time_zone_source = checkFildString(object, "time_zone_source");
        String send_mod_cоnditions = checkFildString(object, "send_mod_cоnditions");
        String send_action_conditions = checkFildString(object, "send_action_conditions");
        String send_radius_conditions = checkFildString(object, "send_radius_conditions");
        Timestamp date_from_condition = checkFildTimestamp(object, "date_from_condition");
        Timestamp date_to_condition = checkFildTimestamp(object, "date_to_condition");
                
        String sender_name = checkFildString(object, "sender_name");
        String message_topic = checkFildString(object, "message_topic");
        String message_channels = checkFildString(object, "message_channels");
        String exception_begin = checkFildString(object, "exception_begin");
        String exception_end = checkFildString(object, "exception_end");
        String message_greeting = checkFildString(object, "message_greeting");
        String message_head = checkFildString(object, "message_head");
        String message_text = checkFildString(object, "message_text");
        
        String action_footer_text = checkFildString(object, "action_footer_text");
        String action_footer_text_link = checkFildString(object, "action_footer_text_link");
        
        String button_link = checkFildString(object, "button_link");
        String button_text = checkFildString(object, "button_text");
        String button_text_link = checkFildString(object, "button_text_link");
        
        String message_warning = checkFildString(object, "message_warning");
        String message_warning_text_link = checkFildString(object, "message_warning_text_link");
        String html = checkFildString(object, "html");
        
        String social_text = checkFildString(object, "social_text");
        String vk_text = checkFildString(object, "vk_text");
        String vk_links = checkFildString(object, "vk_links");
        
        String vk_preview = checkFildString(object, "vk_preview");
        String facebook_text = checkFildString(object, "facebook_text");
        String facebook_links = checkFildString(object, "facebook_links");
        
        String facebook_preview = checkFildString(object, "facebook_preview");
        String telegram_text = checkFildString(object, "telegram_text");
        String telegram_links = checkFildString(object, "telegram_links");
        
        String telegram_preview = checkFildString(object, "telegram_preview");
        String messengers_text = checkFildString(object, "messengers_text");
        String messengers_links = checkFildString(object, "messengers_links");
        
        String messengers_preview = checkFildString(object, "messengers_preview");
        String sms_text = checkFildString(object, "sms_text");
        String sms_links = checkFildString(object, "sms_links");
        
        String sms_preview = checkFildString(object, "sms_preview");
        String count_jump = checkFildString(object, "count_jump");
        String unsubscribe_URL = checkFildString(object, "unsubscribe_URL");
        
        
        ArrayList<ViewSingleMessage> send_condition_mails = new ArrayList<ViewSingleMessage>();
        JsonArray arr0 = checkFildArrayJson(object,"send_condition_mails");
        for(int i = 0; i < arr0.size();i++) {
        send_condition_mails.add(new ViewSingleMessage(checkFildInt((JsonObject)arr0.get(i), "id"), checkFildInt((JsonObject)arr0.get(i), "user_id"), 
                checkFildString((JsonObject)arr0.get(i), "message_head"), checkFildTimestamp((JsonObject)arr0.get(i), "date_send"),
                checkFildTimestamp((JsonObject)arr0.get(i), "date_send").toString().split(" ", 2)[0], 
                checkFildTimestamp((JsonObject)arr0.get(i), "date_send").toString().split(" ", 2)[1], 
                checkFildString((JsonObject)arr0.get(i), "amount_receivers"), checkFildString((JsonObject)arr0.get(i), "amount_opening"),
                checkFildString((JsonObject)arr0.get(i), "amount_visiters"), checkFildString((JsonObject)arr0.get(i),"unknown_field")));
        }
        ArrayList<ViewGroupContacts> sending_groups = new ArrayList<ViewGroupContacts>();
        JsonArray arr1 = checkFildArrayJson(object,"sending_groups");
        for(int i = 0; i < arr1.size();i++) {
            sending_groups.add(new ViewGroupContacts(checkFildInt((JsonObject)arr1.get(i), "id"), 
                    checkFildInt((JsonObject)arr1.get(i), "user_id"), checkFildInt((JsonObject)arr1.get(i), "type_id"), 
                    checkFildString((JsonObject)arr1.get(i), "group_name"), checkFildString((JsonObject)arr1.get(i), "type_name"), 
                    checkFildInt((JsonObject)arr1.get(i), "amount_contacts"), checkFildInt((JsonObject)arr1.get(i), "activated_contacts")));
        }
        
        ArrayList<ViewInvoice> sending_invoices = new ArrayList<ViewInvoice>();
        JsonArray arr2 = checkFildArrayJson(object,"sending_invoices");
        for(int i = 0; i < arr2.size();i++) {
                ArrayList<ShortProduct> products = new ArrayList<ShortProduct>(); 
                JsonArray arr = checkFildArrayJson((JsonObject)arr2.get(i),"products"); 
                for(int j = 0; j < arr.size();i++) { 
                    products.add(new ShortProduct( checkFildInt(((JsonObject)arr.get(j)), "id") , checkFildString(((JsonObject)arr.get(j)), "product_name") )); 
                }

                sending_invoices.add(new ViewInvoice(checkFildInt((JsonObject)arr2.get(i), "id"),
                checkFildInt((JsonObject)arr2.get(i), "user_id"), checkFildString((JsonObject)arr2.get(i), "email"),
                checkFildString((JsonObject)arr2.get(i), "name"), checkFildString((JsonObject)arr2.get(i), "phone"),
                checkFildDouble((JsonObject)arr2.get(i), "amount"), checkFildTimestamp((JsonObject)arr2.get(i), "date_create"),
                checkFildString((JsonObject)arr2.get(i), "payment_status"), products));
        }
        
        ArrayList<ViewGroupContacts> sending_groups_exception = new ArrayList<ViewGroupContacts>();
        JsonArray arr3 = checkFildArrayJson(object,"sending_groups_exception");
        for(int i = 0; i < arr3.size();i++) {
            sending_groups_exception.add(new ViewGroupContacts(checkFildInt((JsonObject)arr3.get(i), "id"), 
                    checkFildInt((JsonObject)arr3.get(i), "user_id"), checkFildInt((JsonObject)arr3.get(i), "type_id"), 
                    checkFildString((JsonObject)arr3.get(i), "group_name"), checkFildString((JsonObject)arr3.get(i), "type_name"), 
                    checkFildInt((JsonObject)arr3.get(i), "amount_contacts"), checkFildInt((JsonObject)arr3.get(i), "activated_contacts")));
        }
        
        ArrayList<ViewInvoice> sending_invoices_exception = new ArrayList<ViewInvoice>();
        JsonArray arr4 = checkFildArrayJson(object,"sending_groups_exception");
        for(int i = 0; i < arr4.size();i++) {
            ArrayList<ShortProduct> products = new ArrayList<ShortProduct>(); 
            JsonArray arr = checkFildArrayJson((JsonObject)arr4.get(i),"products"); 
            for(int j = 0; j < arr.size();i++) { 
                products.add(new ShortProduct( checkFildInt(((JsonObject)arr.get(j)), "id") , checkFildString(((JsonObject)arr.get(j)), "product_name"))); 
                
            }
            
            sending_invoices_exception.add(new ViewInvoice(checkFildInt((JsonObject)arr4.get(i), "id"),
            checkFildInt((JsonObject)arr4.get(i), "user_id"), checkFildString((JsonObject)arr4.get(i), "email"),
            checkFildString((JsonObject)arr4.get(i), "name"), checkFildString((JsonObject)arr4.get(i), "phone"),
            checkFildDouble((JsonObject)arr4.get(i), "amount"), checkFildTimestamp((JsonObject)arr4.get(i), "date_create"),
            checkFildString((JsonObject)arr4.get(i), "payment_status"), products));
        }
        
        ArrayList<MessageLinks> message_links = new ArrayList<MessageLinks>();
        JsonArray arr5 = checkFildArrayJson(object,"message_links");
        for(int i = 0; i < arr5.size(); i++){
            message_links.add(new MessageLinks(checkFildInt((JsonObject)arr5.get(i), "id"),
            checkFildInt((JsonObject)arr5.get(i), "user_id"),
            checkFildInt((JsonObject)arr5.get(i), "message_id"),
            checkFildString((JsonObject)arr5.get(i), "link_address"),
            checkFildString((JsonObject)arr5.get(i), "link_name"),
            checkFildString((JsonObject)arr5.get(i), "link_text_link")));
        }
        
        ArrayList<MessageImageLinks> image_links = new ArrayList<MessageImageLinks>();
        JsonArray arr6 = checkFildArrayJson(object,"image_links");
        for(int i = 0; i < arr6.size();i++) {
            image_links.add(new MessageImageLinks(checkFildInt((JsonObject)arr6.get(i), "id"),
            checkFildInt((JsonObject)arr6.get(i), "user_id"),
            checkFildInt((JsonObject)arr6.get(i), "message_id"),
            checkFildString((JsonObject)arr6.get(i), "image_address"),
            checkFildString((JsonObject)arr6.get(i), "image_name"),
            checkFildString((JsonObject)arr6.get(i), "image_text_link")));     
        }
        
        ArrayList<ViewGroupContacts> message = new ArrayList<ViewGroupContacts>();
        JsonArray arr7 = checkFildArrayJson(object,"mesage");
        for(int i = 0; i < arr7.size();i++) {
            message.add(new ViewGroupContacts(checkFildInt((JsonObject)arr7.get(i), "id"), 
                    checkFildInt((JsonObject)arr7.get(i), "user_id"), checkFildInt((JsonObject)arr7.get(i), "type_id"), 
                    checkFildString((JsonObject)arr7.get(i), "group_name"), checkFildString((JsonObject)arr7.get(i), "type_name"), 
                    checkFildInt((JsonObject)arr7.get(i), "amount_contacts"), checkFildInt((JsonObject)arr7.get(i), "activated_contacts")));
        }
        
        return new SingleMessagesModel(id, user_id, status, date_send, view_date, view_time,
        time_zone, time_zone_source, send_mod_cоnditions, 
        send_action_conditions, send_radius_conditions,
        date_from_condition, date_to_condition, sender_name,
        message_topic, message_channels, exception_begin, exception_end, message_greeting,
        message_head, message_text, action_footer_text,
        action_footer_text_link,button_link, button_text, button_text_link,
        message_warning, message_warning_text_link,
        html, social_text, vk_text, vk_links, vk_preview, facebook_text,
        facebook_links, facebook_preview, telegram_text, telegram_links,
        telegram_preview, messengers_text, messengers_links, messengers_preview,
        sms_text, sms_links, sms_preview, count_jump, unsubscribe_URL,
        send_condition_mails, sending_groups, sending_invoices,
        sending_groups_exception, sending_invoices_exception,
        message_links,image_links, message );
    }
}