package SQLClass.AnswerGroupListObject;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;

public class AnswerGroupListObject {
    int count;
    JsonArray json_message;
    
    public AnswerGroupListObject(int count, JsonArray json_message){
        this.json_message = json_message;
        this.count = count;
    }

    public int getCount() {
        return count;
    }

    public void setCount(int count) {
        this.count = count;
    }

    public JsonArray getJson_message() {
        return json_message;
    }

    public void setJson_message(JsonArray json_message) {
        this.json_message = json_message;
    }
    
    
}
