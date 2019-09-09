package SQLClass.Message;

import MyInterface.ModelInterface;
import SQLClass.Contact.ViewGroupContacts;
import java.sql.Timestamp;
import java.util.ArrayList;

public class SingleMessagesModel implements ModelInterface{
    int id = -1; 
    int user_id = -1;
    String status = "";

    Timestamp date_send; 
    String view_date = "";
    String view_time = "";
    String time_zone = "";
    String time_zone_source = "";

    String send_mod_cоnditions = "";

    String send_action_conditions = "";
    String send_radius_conditions = "";
    Timestamp date_from_condition;
    Timestamp date_to_condition;

    String sender_name = "";
    String message_topic = "";
    String message_channels = "";

    String exception_begin = ""; 
    String exception_end = "";

    String message_greeting = "";
    String message_head = "";
    String message_text = "";

    String action_footer_text = "";
    String action_footer_text_link = "";

    String button_link = "";
    String button_text = "";
    String button_text_link = "";
    
    String message_warning = "";
    String message_warning_text_link = "";

    String html = "";

    String social_text = "";

    String vk_text = "";
    String vk_links = "";
    String vk_preview = "";

    String facebook_text = "";
    String facebook_links = "";
    String facebook_preview = "";

    String telegram_text = "";
    String telegram_links = "";
    String telegram_preview = "";

    String messengers_text = ""; 
    String messengers_links = "";
    String messengers_preview = "";

    String sms_text = "";
    String sms_links = "";
    String sms_preview = "";

    String count_jump = "";
    String unsubscribe_URL = "";

    ArrayList<ViewSingleMessage> send_condition_mails; 
    ArrayList<ViewGroupContacts> sending_groups; 
    ArrayList<ViewInvoice> sending_invoices; 
    ArrayList<ViewGroupContacts> sending_groups_exception; 
    ArrayList<ViewInvoice> sending_invoices_exception;
    ArrayList<MessageLinks> message_links; 
    ArrayList<MessageImageLinks> image_links;
    ArrayList<ViewGroupContacts> message; 

    public SingleMessagesModel(int id, int user_id, String status, Timestamp date_send, String view_date, String view_time,
        String time_zone, String time_zone_source, String send_mod_cоnditions, 
        String send_action_conditions, String send_radius_conditions,
        Timestamp date_from_condition, Timestamp date_to_condition, String sender_name,
        String message_topic, String message_channels, String exception_begin, String exception_end, String message_greeting,
        String message_head, String message_text, String action_footer_text,
        String action_footer_text_link, String button_link, String button_text,
        String button_text_link, String message_warning, String message_warning_text_link,
        String html, String social_text, String vk_text, String vk_links, String vk_preview, String facebook_text,
        String facebook_links, String facebook_preview, String telegram_text, String telegram_links,
        String telegram_preview, String messengers_text, String messengers_links, String messengers_preview,
        String sms_text, String sms_links, String sms_preview, String count_jump, String unsubscribe_URL,
        ArrayList<ViewSingleMessage> send_condition_mails,
        ArrayList<ViewGroupContacts> sending_groups, ArrayList<ViewInvoice> sending_invoices,
        ArrayList<ViewGroupContacts> sending_groups_exception, ArrayList<ViewInvoice> sending_invoices_exception,
        ArrayList<MessageLinks> message_links, ArrayList<MessageImageLinks> image_links, ArrayList<ViewGroupContacts> message)  {

        this.id = id;
        this.user_id = user_id;
        this.status = status;

        this.date_send = date_send; 
        this.view_date = view_date;
        this.view_time = view_time;
        this.time_zone = time_zone;
        this.time_zone_source = time_zone_source;

        this.send_mod_cоnditions = send_mod_cоnditions;

        this.send_action_conditions = send_action_conditions;
        this.send_radius_conditions = send_radius_conditions;
        this.date_from_condition = date_from_condition;
        this.date_to_condition = date_to_condition;

        this.sender_name = sender_name;
        this.message_topic = message_topic;
        this.message_channels = message_channels;

        this.exception_begin = exception_begin; 
        this.exception_end = exception_end;

        this.message_greeting = message_greeting;
        
        this.message_head = message_head;
        this.message_text = message_text;

        this.action_footer_text = action_footer_text;
        this.action_footer_text_link = action_footer_text_link;
        
        this.button_link = button_link;
        this.button_text = button_text;
        this.button_text_link = button_text_link;

        this.message_warning = message_warning;
        this.message_warning_text_link = message_warning_text_link;

        this.html = html;

        this.social_text = social_text;

        this.vk_text = vk_text;
        this.vk_links = vk_links;
        this.vk_preview = vk_preview;

        this.facebook_text = facebook_text;
        this.facebook_links = facebook_links;
        this.facebook_preview = facebook_preview;

        this.telegram_text = telegram_text;
        this.telegram_links = telegram_links;
        this.telegram_preview = telegram_preview;

        this.messengers_text = messengers_text; 
        this.messengers_links = messengers_links;
        this.messengers_preview = messengers_preview;

        this.sms_text = sms_text;
        this.sms_links = sms_links;
        this.sms_preview = sms_preview;

        this.count_jump = count_jump;
        this.unsubscribe_URL = unsubscribe_URL;

        this.send_condition_mails = send_condition_mails; 
        this.sending_groups = sending_groups; 
        this.sending_invoices = sending_invoices; 
        this.sending_groups_exception = sending_groups_exception; 
        this.sending_invoices_exception = sending_invoices_exception;
        this.message_links = message_links; 
        this.image_links = image_links;
        this.message = message;
    }

    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), status, date_send.toString(), view_date,
            view_time, time_zone , time_zone_source, send_mod_cоnditions, send_action_conditions, send_radius_conditions,
            date_from_condition.toString(), date_to_condition.toString(), sender_name, message_topic,
            message_channels, exception_begin, exception_end, message_greeting, message_head, message_text, action_footer_text, action_footer_text_link,
            button_link, button_text, button_text_link,
            message_warning, message_warning_text_link, html, social_text, vk_text, vk_links, vk_preview, facebook_text,
            facebook_links, facebook_preview, telegram_text, telegram_links, telegram_preview, messengers_text,
            messengers_links, messengers_preview, sms_text, sms_links, sms_preview, count_jump, unsubscribe_URL
        };
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id),
            String.valueOf(user_id), 
            "'" + status.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + date_send.toString() + "'", 
            "'" + view_date.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + view_time.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + time_zone.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + time_zone_source.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + send_mod_cоnditions.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + send_action_conditions.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + send_radius_conditions.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + date_from_condition.toString() + "'", 
            "'" + date_to_condition.toString() + "'", 
            "'" + sender_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + message_channels.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + exception_begin.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + exception_end.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + message_greeting.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + message_head.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + message_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + action_footer_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + action_footer_text_link.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + button_link.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + button_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + button_text_link.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",   
            "'" + message_warning.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + message_warning_text_link.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + html.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + social_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + vk_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + vk_links.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + vk_preview.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + facebook_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + facebook_links.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + facebook_preview.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + telegram_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + telegram_links.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + telegram_preview.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + messengers_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + messengers_links.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + messengers_preview.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + sms_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + sms_links.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + sms_preview.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + count_jump.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + unsubscribe_URL.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "status", "date_send", "view_date", "view_time", "time_zone", 
            "time_zone_source", "send_mod_cоnditions", "send_action_conditions", "send_radius_conditions", 
            "date_from_condition", "date_to_condition", "sender_name", 
            "message_channels", "exception_begin", "exception_end", "message_greeting", "message_head", "message_text", 
            "action_footer_text", "action_footer_text_link", "button_link", "button_text", "button_text_link",
            "message_warning", "message_warning_text_link", "html", "social_text", "vk_text", 
            "vk_links", "vk_preview", "facebook_text", "facebook_links", "facebook_preview", "telegram_text",
            "telegram_links", "telegram_preview", 
            "messengers_text", "messengers_links", "messengers_preview", "sms_text", "sms_links", "sms_preview",
            "count_jump", "unsubscribe_URL"}; 
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

    public String getStatus() {
        return status;
    }

    public void setStatus(String status) {
        this.status = status;
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

    public String getTime_zone_source() {
        return time_zone_source;
    }

    public void setTime_zone_source(String time_zone_source) {
        this.time_zone_source = time_zone_source;
    }

    public String getSend_mod_cоnditions() {
        return send_mod_cоnditions;
    }

    public void setSend_mod_cоnditions(String send_mod_cоnditions) {
        this.send_mod_cоnditions = send_mod_cоnditions;
    }

    public String getSend_action_conditions() {
        return send_action_conditions;
    }

    public void setSend_action_conditions(String send_action_conditions) {
        this.send_action_conditions = send_action_conditions;
    }

    public String getSend_radius_conditions() {
        return send_radius_conditions;
    }

    public void setSend_radius_conditions(String send_radius_conditions) {
        this.send_radius_conditions = send_radius_conditions;
    }

    public Timestamp getDate_from_condition() {
        return date_from_condition;
    }

    public void setDate_from_condition(Timestamp date_from_condition) {
        this.date_from_condition = date_from_condition;
    }

    public Timestamp getDate_to_condition() {
        return date_to_condition;
    }

    public void setDate_to_condition(Timestamp date_to_condition) {
        this.date_to_condition = date_to_condition;
    }

    public ArrayList<ViewSingleMessage> getSend_condition_mails() {
        return send_condition_mails;
    }

    public void setSend_condition_mails(ArrayList<ViewSingleMessage>  send_condition_mails) {
        this.send_condition_mails = send_condition_mails;
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

    public String getMessage_channels() {
        return message_channels;
    }

    public void setMessage_channels(String message_channels) {
        this.message_channels = message_channels;
    }

    public String getException_begin() {
        return exception_begin;
    }

    public void setException_begin(String exception_begin) {
        this.exception_begin = exception_begin;
    }

    public String getException_end() {
        return exception_end;
    }

    public void setException_end(String exception_end) {
        this.exception_end = exception_end;
    }

    public String getMessage_head() {
        return message_head;
    }

    public void setMessage_head(String message_head) {
        this.message_head = message_head;
    }

    public String getMessage_text() {
        return message_text;
    }

    public void setMessage_text(String message_text) {
        this.message_text = message_text;
    }

    public String getAction_footer_text() {
        return action_footer_text;
    }

    public void setAction_footer_text(String action_footer_text) {
        this.action_footer_text = action_footer_text;
    }

    public String getAction_footer_text_link() {
        return action_footer_text_link;
    }

    public void setAction_footer_text_link(String action_footer_text_link) {
        this.action_footer_text_link = action_footer_text_link;
    }

    public String getButton_link() {
        return button_link;
    }

    public void setButton_link(String button_link) {
        this.button_link = button_link;
    }

    public String getButton_text() {
        return button_text;
    }

    public void setButton_text(String button_text) {
        this.button_text = button_text;
    }

    public String getButton_text_link() {
        return button_text_link;
    }

    public void setButton_text_link(String button_text_link) {
        this.button_text_link = button_text_link;
    }

    public String getMessage_warning() {
        return message_warning;
    }

    public void setMessage_warning(String message_warning) {
        this.message_warning = message_warning;
    }

    public String getMessage_warning_text_link() {
        return message_warning_text_link;
    }

    public void setMessage_warning_text_link(String message_warning_text_link) {
        this.message_warning_text_link = message_warning_text_link;
    }

    public String getHtml() {
        return html;
    }

    public void setHtml(String html) {
        this.html = html;
    }

    public String getSocial_text() {
        return social_text;
    }

    public void setSocial_text(String social_text) {
        this.social_text = social_text;
    }

    public String getVk_text() {
        return vk_text;
    }

    public void setVk_text(String vk_text) {
        this.vk_text = vk_text;
    }

    public String getVk_links() {
        return vk_links;
    }

    public void setVk_links(String vk_links) {
        this.vk_links = vk_links;
    }

    public String getVk_preview() {
        return vk_preview;
    }

    public void setVk_preview(String vk_preview) {
        this.vk_preview = vk_preview;
    }

    public String getFacebook_text() {
        return facebook_text;
    }

    public void setFacebook_text(String facebook_text) {
        this.facebook_text = facebook_text;
    }

    public String getFacebook_links() {
        return facebook_links;
    }

    public void setFacebook_links(String facebook_links) {
        this.facebook_links = facebook_links;
    }

    public String getFacebook_preview() {
        return facebook_preview;
    }

    public void setFacebook_preview(String facebook_preview) {
        this.facebook_preview = facebook_preview;
    }

    public String getTelegram_text() {
        return telegram_text;
    }

    public void setTelegram_text(String telegram_text) {
        this.telegram_text = telegram_text;
    }

    public String getTelegram_links() {
        return telegram_links;
    }

    public void setTelegram_links(String telegram_links) {
        this.telegram_links = telegram_links;
    }

    public String getTelegram_preview() {
        return telegram_preview;
    }

    public void setTelegram_preview(String telegram_preview) {
        this.telegram_preview = telegram_preview;
    }

    public String getMessengers_text() {
        return messengers_text;
    }

    public void setMessengers_text(String messengers_text) {
        this.messengers_text = messengers_text;
    }

    public String getMessengers_links() {
        return messengers_links;
    }

    public void setMessengers_links(String messengers_links) {
        this.messengers_links = messengers_links;
    }

    public String getMessengers_preview() {
        return messengers_preview;
    }

    public void setMessengers_preview(String messengers_preview) {
        this.messengers_preview = messengers_preview;
    }

    public String getSms_text() {
        return sms_text;
    }

    public void setSms_text(String sms_text) {
        this.sms_text = sms_text;
    }

    public String getSms_links() {
        return sms_links;
    }

    public void setSms_links(String sms_links) {
        this.sms_links = sms_links;
    }

    public String getSms_preview() {
        return sms_preview;
    }

    public void setSms_preview(String sms_preview) {
        this.sms_preview = sms_preview;
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

    public ArrayList<ViewGroupContacts> getSending_groups() {
        return sending_groups;
    }

    public void setSending_groups(ArrayList<ViewGroupContacts> sending_groups) {
        this.sending_groups = sending_groups;
    }

    public ArrayList<ViewInvoice> getSending_invoices() {
        return sending_invoices;
    }

    public void setSending_invoices(ArrayList<ViewInvoice> sending_invoices) {
        this.sending_invoices = sending_invoices;
    }

    public ArrayList<ViewGroupContacts> getSending_groups_exception() {
        return sending_groups_exception;
    }

    public void setSending_groups_exception(ArrayList<ViewGroupContacts> sending_groups_exception) {
        this.sending_groups_exception = sending_groups_exception;
    }

    public ArrayList<ViewInvoice> getSending_invoices_exception() {
        return sending_invoices_exception;
    }

    public void setSending_invoices_exception(ArrayList<ViewInvoice> sending_invoices_exception) {
        this.sending_invoices_exception = sending_invoices_exception;
    }

    public ArrayList<MessageLinks> getMessage_links() {
        return message_links;
    }

    public void setMessage_links(ArrayList<MessageLinks> message_links) {
        this.message_links = message_links;
    }

    public ArrayList<MessageImageLinks> getImage_links() {
        return image_links;
    }

    public void setImage_links(ArrayList<MessageImageLinks> image_links) {
        this.image_links = image_links;
    }

    public ArrayList<ViewGroupContacts> getMessage() {
        return message;
    }

    public void setMessage(ArrayList<ViewGroupContacts> message) {
        this.message = message;
    }
}
