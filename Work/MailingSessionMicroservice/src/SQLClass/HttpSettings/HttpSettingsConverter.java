package SQLClass.HttpSettings;

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
import java.util.ArrayList;


public class HttpSettingsConverter implements JsonSerializer<HttpSettings>, JsonDeserializer<HttpSettings>, InterfaceConverter{
    @Override
    public JsonElement serialize(HttpSettings src, Type type, JsonSerializationContext jsc) {
        JsonObject settings = new JsonObject();
        settings.addProperty("id", src.getId());
        settings.addProperty("user_id", src.getUser_id());
        
        settings.addProperty("emails", src.getEmails());
        settings.addProperty("url_sub", src.getUrl_sub());
        settings.addProperty("use_default_act", src.getUse_default_act());
        settings.addProperty("url_act", src.getUrl_act());
        settings.addProperty("use_default_whitelist", src.getUse_default_whitelist());
        settings.addProperty("url_unsub", src.getUrl_unsub());
        settings.addProperty("use_default_panel", src.getUse_default_panel());
        settings.addProperty("auto_unsub", src.getAuto_unsub());
        
        settings.addProperty("logo", src.getLogo());
        settings.addProperty("logo_link", src.getLogo_link());
        settings.addProperty("description", src.getDescription());
        settings.addProperty("contact_icon", src.getContact_icon());
        settings.addProperty("main_contact", src.getMain_contact());
        settings.addProperty("main_contact_link", src.getMain_contact_link());
        settings.addProperty("additional_contact", src.getAdditional_contact());
        settings.addProperty("additional_contact_link", src.getAdditional_contact_link());
        settings.addProperty("picture", src.getPicture());
        settings.addProperty("greeting", src.getGreeting());
        settings.addProperty("template", src.getTemplate());
        settings.addProperty("social_message", src.getSocial_message());
        JsonArray icons = new JsonArray();
        for(int i = 0; i < src.getSocial_icons().size();i++) {
            JsonObject icon = new JsonObject();
            icon.addProperty("id", src.getSocial_icons().get(i).getId());
            icon.addProperty("user_id", src.getSocial_icons().get(i).getUser_id());
            icon.addProperty("http_settings_id", src.getSocial_icons().get(i).getHttp_settings_id());
            icon.addProperty("icon_path", src.getSocial_icons().get(i).getIcon_path());
            icon.addProperty("link", src.getSocial_icons().get(i).getLink());
            icons.add(icon);
        }
        settings.add("social_icons", icons);
        JsonArray user_emails = new JsonArray();
        for(int i = 0; i < src.getUser_emails().size();i++) {
            JsonObject user_email = new JsonObject();
            user_email.addProperty("id", src.getUser_emails().get(i).getId());
            user_email.addProperty("user_id", src.getUser_emails().get(i).getUser_id());
            user_email.addProperty("http_settings_id", src.getUser_emails().get(i).getHttp_settings_id());
            user_email.addProperty("sender_name", src.getUser_emails().get(i).getSender_name());
            user_email.addProperty("sender_email", src.getUser_emails().get(i).getSender_email());
            user_email.addProperty("mailing_name", src.getUser_emails().get(i).getMailing_name());
            user_email.addProperty("main_url", src.getUser_emails().get(i).getMain_url());
            user_emails.add(user_email);
        }
        settings.add("user_emails", user_emails);
        return settings;
    }

    @Override
    public HttpSettings deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id"); 
        int user_id = checkFildInt(object, "user_id");
        
        String emails = checkFildString(object, "emails");
        String url_sub = checkFildString(object, "url_sub");
        String use_default_act = checkFildString(object, "use_default_act");
        String url_act = checkFildString(object, "url_act");
        String use_default_whitelist = checkFildString(object, "use_default_whitelist");
        String url_unsub = checkFildString(object, "url_unsub");
        String use_default_panel = checkFildString(object, "use_default_panel");
        String auto_unsub = checkFildString(object, "auto_unsub");
        
        String logo = checkFildString(object, "logo");
        String logo_link = checkFildString(object, "logo_link");
        String description = checkFildString(object, "description");
        String contact_icon = checkFildString(object, "contact_icon");
        String main_contact = checkFildString(object, "main_contact");
        String main_contact_link = checkFildString(object, "main_contact_link");
        String additional_contact = checkFildString(object, "additional_contact");
        String additional_contact_link = checkFildString(object, "additional_contact_link");
        String picture = checkFildString(object, "picture");
        String greeting = checkFildString(object, "greeting");
        String template = checkFildString(object, "template");
        String social_message = checkFildString(object, "social_message");

        ArrayList<SocialIcons> socialIcons = new ArrayList<SocialIcons>();
        JsonArray arr1 = checkFildArrayJson(object,"social_icons");
        for(int i = 0; i < arr1.size();i++) {
            socialIcons.add(new SocialIcons(checkFildInt(((JsonObject)arr1.get(i)), "id"),
                    checkFildInt(((JsonObject)arr1.get(i)), "user_id"), 
                    checkFildInt(((JsonObject)arr1.get(i)), "http_settings_id"), 
                    checkFildString(((JsonObject)arr1.get(i)), "icon_path"),
                    checkFildString(((JsonObject)arr1.get(i)), "link")));
        }
        
        ArrayList<UserEmails> user_emails = new ArrayList<UserEmails>();
        JsonArray arr2 = checkFildArrayJson(object,"user_emails");
        for(int i = 0; i < arr2.size();i++) {
            user_emails.add(new UserEmails(checkFildInt(((JsonObject)arr2.get(i)), "id"),
                    checkFildInt(((JsonObject)arr2.get(i)), "user_id"),
                    checkFildInt(((JsonObject)arr2.get(i)), "http_settings_id"),
                    checkFildString(((JsonObject)arr2.get(i)), "sender_name"),
                    checkFildString(((JsonObject)arr2.get(i)), "sender_email"),
                    checkFildString(((JsonObject)arr2.get(i)), "mailing_name"),
                    checkFildString(((JsonObject)arr2.get(i)), "main_url")));
        }
        

        return new HttpSettings(id, user_id, emails, url_sub, use_default_act, url_act, use_default_whitelist,
                url_unsub, use_default_panel, auto_unsub, logo, logo_link, description, 
            contact_icon, main_contact, main_contact_link, additional_contact, additional_contact_link, picture,
            greeting, template, social_message, socialIcons, user_emails);
    }
}
