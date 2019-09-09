package SQLClass.HttpSettings;

import MyInterface.InterfaceConverter;
import com.google.gson.JsonDeserializationContext;
import com.google.gson.JsonDeserializer;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;
import java.lang.reflect.Type;

public class UserEmailsConverter implements JsonSerializer<UserEmails>, JsonDeserializer<UserEmails>, InterfaceConverter{

    @Override
    public JsonElement serialize(UserEmails src, Type type, JsonSerializationContext jsc) {
        JsonObject mails = new JsonObject();
        mails.addProperty("id", src.getId());
        mails.addProperty("user_id", src.getUser_id());
        mails.addProperty("http_settings_id", src.getHttp_settings_id());
        mails.addProperty("sender_name", src.getSender_name());
        mails.addProperty("sender_email", src.getSender_email());
        mails.addProperty("mailing_name", src.getMailing_name());
        mails.addProperty("main_url", src.getMain_url());
        return mails;
    }

    @Override
    public UserEmails deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        int http_settings_id = checkFildInt(object, "http_settings_id");
        String sender_name = checkFildString(object, "sender_name");
        String sender_email = checkFildString(object, "sender_email");
        String mailing_name = checkFildString(object, "mailing_name");
        String main_url = checkFildString(object, "main_url");
        return new UserEmails(id,user_id, http_settings_id, sender_name, sender_email, mailing_name, main_url);
    }
}
