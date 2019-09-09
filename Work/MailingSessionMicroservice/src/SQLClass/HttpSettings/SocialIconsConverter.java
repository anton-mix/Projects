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

public class SocialIconsConverter implements JsonSerializer<SocialIcons>, JsonDeserializer<SocialIcons>, InterfaceConverter{

    @Override
    public JsonElement serialize(SocialIcons src, Type type, JsonSerializationContext jsc) {
        JsonObject icons = new JsonObject();
        icons.addProperty("id", src.getId());
        icons.addProperty("user_id", src.getUser_id());
        icons.addProperty("http_settings_id", src.getHttp_settings_id());
        icons.addProperty("icon_path", src.getIcon_path());
        icons.addProperty("link", src.getLink());
        return icons;
    }

    @Override
    public SocialIcons deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        int id = checkFildInt(object, "id");
        int user_id = checkFildInt(object, "user_id");
        int http_settings_id = checkFildInt(object, "http_settings_id");
        String icon_path = checkFildString(object, "icon_path");
        String link = checkFildString(object, "link");
        return new SocialIcons(id,user_id, http_settings_id, icon_path, link);
    }
    
}
