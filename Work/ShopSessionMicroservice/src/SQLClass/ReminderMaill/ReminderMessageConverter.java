package SQLClass.ReminderMaill;

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
import java.sql.Timestamp;
import java.util.ArrayList;

public class ReminderMessageConverter implements JsonSerializer<ReminderMessage>, JsonDeserializer<ReminderMessage>, InterfaceConverter{
    @Override
    public JsonElement serialize(ReminderMessage src, Type type, JsonSerializationContext jsc) {
        JsonObject message = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            message.addProperty(key[i], values[i]);
        
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
                
        return message;
    }
    
    @Override
    public ReminderMessage deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id"); 
        int reminder_id = checkFildInt(object, "reminder_id"); 
        int user_id = checkFildInt(object, "user_id");
        Timestamp date_send = checkFildTimestamp(object, "date_send");
        String view_date = date_send.toString().split(" ", 2)[0];
        String view_time = date_send.toString().split(" ", 2)[1];
        String time_zone = checkFildString(object, "time_zone");
        String time_zone_source = checkFildString(object, "time_zone_source");
                
        String sender_name = checkFildString(object, "sender_name");
        String message_topic = checkFildString(object, "message_topic");
        String message_channels = checkFildString(object, "message_channels");
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
        
        return new ReminderMessage(id, user_id, reminder_id, date_send, view_date, view_time,
        time_zone, time_zone_source, sender_name, message_topic, message_channels, message_greeting,
        message_head, message_text, action_footer_text, action_footer_text_link,button_link, 
        button_text, button_text_link, message_warning, message_warning_text_link,
        html, social_text, vk_text, vk_links, vk_preview, facebook_text,
        facebook_links, facebook_preview, telegram_text, telegram_links,
        telegram_preview, messengers_text, messengers_links, messengers_preview,
        sms_text, sms_links, sms_preview, count_jump, unsubscribe_URL,
        message_links,image_links);
    }
}
