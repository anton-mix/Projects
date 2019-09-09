package ImageObjects;

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

public class LoadImageConverter implements JsonSerializer<LoadImage>, JsonDeserializer<LoadImage>, InterfaceConverter  {
    @Override
    public JsonElement serialize(LoadImage src, Type type, JsonSerializationContext jsc) {
        JsonObject session_object = new JsonObject();
        session_object.addProperty("user_id", src.getUser_id());
        session_object.addProperty("microservice_id", src.getMicroservice_id());
        session_object.addProperty("image_id", src.getImage_id());
        return session_object;
    }
    
    @Override
    public LoadImage deserialize(JsonElement json, Type type, JsonDeserializationContext jdc) throws JsonParseException {
        JsonObject object = json.getAsJsonObject();
        String user_id = checkFildString(object, "user_id");
        int microservice_id = checkFildInt(object, "microservice_id");
        int image_id = checkFildInt(object, "image_id");
        ArrayList<String> part_path = new ArrayList<String>();
        JsonArray part_path_arr = checkFildArrayJson(object, "part_path");
        for(int i = 0; i < part_path_arr.size(); i ++) {
            part_path.add(part_path_arr.get(i).getAsString());
        }
        return new LoadImage(user_id, microservice_id, image_id, part_path);
    }
}
