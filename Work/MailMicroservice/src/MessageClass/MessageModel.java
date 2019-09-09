package MessageClass;

import MyInterface.ModelInterface;
import java.sql.Timestamp;
import java.util.ArrayList;

public class MessageModel implements ModelInterface{
    int id = -1; 
    int user_id = -1;
    Timestamp date_send;
    String view_date = "";
    String view_time = "";
    String time_zone = "";
    String sender_name = "";
    String message_topic = "";
    String message_head = "";
    int message_format = -1;
    String message_text = "";
    String action_header_link = "";
    String action_header_text = ""; 
    String action_footer_link = ""; 
    String action_footer_text = ""; 
    String message_warning = ""; 
    int message_channels = -1;
    Timestamp exception_begin;
    Timestamp exception_end; 
    String count_jump = ""; 
    String unsubscribe_URL = "";
    ArrayList<MessageLinks> message_links;
    ArrayList<ImageLinks> image_links;

    public MessageModel(int id, int user_id, Timestamp date_send, String view_date, String view_time, String time_zone,
    String sender_name, String message_topic, String message_head, int message_format, String message_text,
    String action_header_link, String action_header_text, String action_footer_link,
    String action_footer_text, String message_warning, int message_channels, Timestamp exception_begin,
    Timestamp exception_end, String count_jump, String unsubscribe_URL, ArrayList<MessageLinks> message_links, ArrayList<ImageLinks> image_links) {
        this.id = id;
        this.user_id = user_id;
        this.date_send = date_send;
        this.view_date = view_date;
        this.view_time = view_time;
        this.time_zone = time_zone;
        this.sender_name = sender_name;
        this.message_topic = message_topic;
        this.message_head = message_head;
        this.message_format = message_format;
        this.message_text = message_text;
        this.action_header_link = action_header_link;
        this.action_header_text = action_header_text;
        this.action_footer_link = action_footer_link;
        this.action_footer_text = action_footer_text;
        this.message_warning = message_warning;
        this.message_channels = message_channels;
        this.exception_begin = exception_begin;
        this.exception_end = exception_end;
        this.count_jump = count_jump;
        this.unsubscribe_URL = unsubscribe_URL;
        this.message_links = message_links;
        this.image_links = image_links;
    }
    
    public MessageModel(int user_id, Timestamp date_send, String view_date, String view_time, String time_zone,
    String sender_name, String message_topic, String message_head, int message_format, String message_text,
    String action_header_link, String action_header_text, String action_footer_link,
    String action_footer_text, String message_warning, int message_channels, Timestamp exception_begin,
    Timestamp exception_end, String count_jump, String unsubscribe_URL) {
        this.user_id = user_id;
        this.date_send = date_send;
        this.view_date = view_date;
        this.view_time = view_time;
        this.time_zone = time_zone;
        this.sender_name = sender_name;
        this.message_topic = message_topic;
        this.message_head = message_head;
        this.message_format = message_format;
        this.message_text = message_text;
        this.action_header_link = action_header_link;
        this.action_header_text = action_header_text;
        this.action_footer_link = action_footer_link;
        this.action_footer_text = action_footer_text;
        this.message_warning = message_warning;
        this.message_channels = message_channels;
        this.exception_begin = exception_begin;
        this.exception_end = exception_end;
        this.count_jump = count_jump;
        this.unsubscribe_URL = unsubscribe_URL;
        this.message_links = new ArrayList<MessageLinks>();
        this.image_links = new ArrayList<ImageLinks>();
    }
        
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), date_send.toString(), view_date, view_time, 
            time_zone, sender_name, message_topic, message_head, String.valueOf(message_format), message_text,
            action_header_link, action_header_text, action_footer_link,
            action_footer_text, message_warning, String.valueOf(message_channels), exception_begin.toString(),
            exception_end.toString(), count_jump, unsubscribe_URL};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id),
            "'" + date_send.toString().replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + view_date.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + view_time.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + time_zone.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + sender_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + message_topic.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + message_head.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            String.valueOf(message_format), 
            "'" + message_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + action_header_link.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + action_header_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + action_footer_link.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + action_footer_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + message_warning.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            String.valueOf(message_channels),
            "'" + exception_begin.toString().replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + exception_end.toString().replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + count_jump.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + unsubscribe_URL.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "date_send", "view_date", "view_time", 
            "time_zone", "sender_name", "message_topic", "message_head", "message_format", "message_text",
            "action_header_link", "action_header_text", "action_footer_link",
            "action_footer_text", "message_warning", "message_channels", "exception_begin",
            "exception_end", "count_jump", "unsubscribe_URL"}; 
    }

    @Override
    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public Timestamp getDate_send() {
        return date_send;
    }

    public void setDate_send(Timestamp date_send) {
        this.date_send = date_send;
    }

    public String getView_date() {
        return view_date;
    }

    public void setView_date(String view_date) {
        this.view_date = view_date;
    }

    public String getView_time() {
        return view_time;
    }

    public void setView_time(String view_time) {
        this.view_time = view_time;
    }

    public String getTime_zone() {
        return time_zone;
    }

    public void setTime_zone(String time_zone) {
        this.time_zone = time_zone;
    }

    public String getSender_name() {
        return sender_name;
    }

    public void setSender_name(String sender_name) {
        this.sender_name = sender_name;
    }

    public String getMessage_topic() {
        return message_topic;
    }

    public void setMessage_topic(String message_topic) {
        this.message_topic = message_topic;
    }

    public String getMessage_head() {
        return message_head;
    }

    public void setMessage_head(String message_head) {
        this.message_head = message_head;
    }

    public int getMessage_format() {
        return message_format;
    }

    public void setMessage_format(int message_format) {
        this.message_format = message_format;
    }

    public String getMessage_text() {
        return message_text;
    }

    public void setMessage_text(String message_text) {
        this.message_text = message_text;
    }

    public String getAction_header_link() {
        return action_header_link;
    }

    public void setAction_header_link(String action_header_link) {
        this.action_header_link = action_header_link;
    }

    public String getAction_header_text() {
        return action_header_text;
    }

    public void setAction_header_text(String action_header_text) {
        this.action_header_text = action_header_text;
    }

    public String getAction_footer_link() {
        return action_footer_link;
    }

    public void setAction_footer_link(String action_footer_link) {
        this.action_footer_link = action_footer_link;
    }

    public String getAction_footer_text() {
        return action_footer_text;
    }

    public void setAction_footer_text(String action_footer_text) {
        this.action_footer_text = action_footer_text;
    }

    public String getMessage_warning() {
        return message_warning;
    }

    public void setMessage_warning(String message_warning) {
        this.message_warning = message_warning;
    }

    public int getMessage_channels() {
        return message_channels;
    }

    public void setMessage_channels(int message_channels) {
        this.message_channels = message_channels;
    }

    public Timestamp getException_begin() {
        return exception_begin;
    }

    public void setException_begin(Timestamp exception_begin) {
        this.exception_begin = exception_begin;
    }

    public Timestamp getException_end() {
        return exception_end;
    }

    public void setException_end(Timestamp exception_end) {
        this.exception_end = exception_end;
    }

    public String getCount_jump() {
        return count_jump;
    }

    public void setCount_jump(String count_jump) {
        this.count_jump = count_jump;
    }

    public String getUnsubscribe_URL() {
        return unsubscribe_URL;
    }

    public void setUnsubscribe_URL(String unsubscribe_URL) {
        this.unsubscribe_URL = unsubscribe_URL;
    }

    public ArrayList<MessageLinks> getMessage_links() {
        return message_links;
    }

    public void setMessage_links(ArrayList<MessageLinks> message_links) {
        this.message_links = message_links;
    }

    public ArrayList<ImageLinks> getImage_links() {
        return image_links;
    }

    public void setImage_links(ArrayList<ImageLinks> image_links) {
        this.image_links = image_links;
    }
}
