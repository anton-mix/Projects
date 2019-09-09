
package MessageClass;

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


public class MessageModelConverter implements JsonSerializer<MessageModel>, JsonDeserializer<MessageModel>, InterfaceConverter{
    @Override
    public JsonElement serialize(MessageModel src, Type type, JsonSerializationContext jsc) {
        JsonObject message = new JsonObject();
        String[] key = src.getKeyArray();
        String[] values = src.getArray();
        for(int i = 0; i < values.length; i++) 
            message.addProperty(key[i], values[i]);
                
        JsonArray message_links = new JsonArray();
        for(int i = 0; i < src.message_links.size();i++) {
            JsonObject message_link = new JsonObject();
            message_link.addProperty("link_address", src.message_links.get(i).getLink_address());
            message_link.addProperty("link_name", src.message_links.get(i).getLink_name());
            message_links.add(message_link);
        }
        message.add("message_links", message_links);
        
        JsonArray image_links = new JsonArray();
        for(int i = 0; i < src.image_links.size();i++) {
            JsonObject image_link = new JsonObject();
            image_link.addProperty("image_address", src.image_links.get(i).getImage_address());
            image_link.addProperty("image_name", src.image_links.get(i).getImage_name());
            image_links.add(image_link);
        }
        message.add("image_links", image_links);
       
        return message;
    }

    @Override
    public MessageModel deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id"); 
        int user_id = checkFildInt(object, "user_id");
        Timestamp date_send = checkFildTimestamp(object, "date_send");
        String view_date = date_send.toString().split(" ", 2)[0];
        String view_time = date_send.toString().split(" ", 2)[1];
        String time_zone = checkFildString(object, "time_zone");
        String sender_name = checkFildString(object, "sender_name");
        String message_topic = checkFildString(object, "message_topic");
        String message_head = checkFildString(object, "message_head");
        int message_format = checkFildInt(object, "message_format");
        String message_text = checkFildString(object, "message_text");
        String action_header_link = checkFildString(object, "action_header_link");
        String action_header_text = checkFildString(object, "action_header_text");
        String action_footer_link = checkFildString(object, "action_footer_link");
        String action_footer_text = checkFildString(object, "action_footer_text");
        String message_warning = checkFildString(object, "message_warning");
        int message_channels = checkFildInt(object, "message_channels");
        Timestamp exception_begin = checkFildTimestamp(object, "exception_begin");
        Timestamp exception_end = checkFildTimestamp(object, "exception_end"); 
        String count_jump = checkFildString(object, "count_jump");
        String unsubscribe_URL = checkFildString(object, "unsubscribe_URL");
        
        ArrayList<MessageLinks> message_links = new ArrayList<MessageLinks>();
        JsonArray arr1 = checkFildArrayJson(object,"message_links");
        for(int i = 0; i < arr1.size();i++) {
            message_links.add(new MessageLinks(((JsonObject)arr1.get(i)).get("link_address").getAsString(), ((JsonObject)arr1.get(i)).get("link_name").getAsString()));
        }
        
        ArrayList<ImageLinks> image_links = new ArrayList<ImageLinks>();
        JsonArray arr2 = checkFildArrayJson(object,"image_links");
        for(int i = 0; i < arr2.size();i++) {
            image_links.add(new ImageLinks(((JsonObject)arr2.get(i)).get("image_address").getAsString(), ((JsonObject)arr2.get(i)).get("image_name").getAsString()));
        }
        
        return new MessageModel(id, user_id, date_send, view_date, view_time, 
            time_zone, sender_name, message_topic, message_head, message_format, message_text,
            action_header_link, action_header_text, action_footer_link,
            action_footer_text, message_warning, message_channels, exception_begin,
            exception_end, count_jump, unsubscribe_URL, message_links, image_links);
    }
}
