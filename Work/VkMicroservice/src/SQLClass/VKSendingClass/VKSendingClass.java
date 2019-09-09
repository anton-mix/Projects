package SQLClass.VKSendingClass;

import java.util.ArrayList;


public class VKSendingClass {
    int owner_id;
    String group_id;
    String user_token;
    String group_token;
    String message_text;
    ArrayList<Integer> user_ids;
    
    public VKSendingClass(){
         
    }
    
    public VKSendingClass(String user_token, String group_token, String message_text, ArrayList<Integer> user_ids, int owner_id, String group_id){
        this.user_ids = user_ids;
        this.user_token = user_token;
        this.message_text = message_text;
        this.owner_id = owner_id;
        this.group_id = group_id;
        this.group_token = group_token;
        
    }

    public String getUser_token() {
        return user_token;
    }

    public void setUser_token(String user_token) {
        this.user_token = user_token;
    }

    public String getGroup_token() {
        return group_token;
    }

    public void setGroup_token(String group_token) {
        this.group_token = group_token;
    }
    
    public String getMessage_text() {
        return message_text;
    }

    public void setMessage_text(String message_text) {
        this.message_text = message_text;
    }

    public ArrayList<Integer> getUser_ids() {
        return user_ids;
    }

    public void setUser_ids(ArrayList<Integer> user_ids) {
        this.user_ids = user_ids;
    }

    public int getOwner_id() {
        return owner_id;
    }

    public void setOwner_id(int owner_id) {
        this.owner_id = owner_id;
    }

    public String getGroup_id() {
        return group_id;
    }

    public void setGroup_id(String group_id) {
        this.group_id = group_id;
    }
}
