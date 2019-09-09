package SQLClass.Contact;

import MyInterface.ModelInterface;
import SQLClass.User.ShortUserEmail;
import java.util.ArrayList;

public class GroupContacts implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int type_id = -1;
    String type_name = "";
    String group_name = ""; 
    int amount_contacts = -1;
    int activated_contacts = -1;
    String activation_text = "";
    String unsubscribe_url = "";
    int sender_id = -1;
    ArrayList<ShortUserEmail> sender;
    
    public GroupContacts(int id, int user_id, String group_name, int type_id, String type_name, int amount_contacts, int activated_contacts,
            String activation_text, String unsubscribe_url, int sender_id){
        this.id = id;
        this.user_id = user_id;
        this.group_name = group_name;
        this.type_id = type_id;
        this.type_name = type_name;
        this.amount_contacts = amount_contacts;
        this.activated_contacts = activated_contacts;
        this.activation_text = activation_text;
        this.unsubscribe_url = unsubscribe_url;
        this.sender_id = sender_id;
        this.sender = new ArrayList<ShortUserEmail>();
    }
    
    public GroupContacts(int user_id, int type_id, String group_name, String type_name, int amount_contacts, int activated_contacts,
            String activation_text, String unsubscribe_url, int sender_id){
        this.user_id = user_id;
        this.type_id = type_id;
        this.group_name = group_name;
        this.type_name = type_name;
        this.amount_contacts = amount_contacts;
        this.activated_contacts = activated_contacts;
        this.activation_text = activation_text;
        this.unsubscribe_url = unsubscribe_url;
        this.sender_id = sender_id;
        this.sender = new ArrayList<ShortUserEmail>();
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(type_id), group_name, String.valueOf(amount_contacts), 
            String.valueOf(activated_contacts), activation_text, unsubscribe_url, String.valueOf(sender_id), type_name};
    }
    
    @Override
    public String[] getSqlArray() {
            return new String[] {
                String.valueOf(id),
                String.valueOf(user_id), 
                String.valueOf(type_id), 
                "'" + group_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
                String.valueOf(amount_contacts), 
                String.valueOf(activated_contacts),
                "'" + activation_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
                "'" + unsubscribe_url.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
                String.valueOf(sender_id)
        };
    }
    
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "type_id", "group_name", "amount_contacts", "activated_contacts", "activation_text",
        "unsubscribe_url", "sender_id", "type_name"}; 
    }

    @Override
    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public int getType_id() {
        return type_id;
    }

    public void setType_id(int type_id) {
        this.type_id = type_id;
    }

    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public String getGroup_name() {
        return group_name;
    }

    public void setGroup_name(String group_name) {
        this.group_name = group_name;
    }

    public int getAmount_contacts() {
        return amount_contacts;
    }

    public void setAmount_contacts(int amount_contacts) {
        this.amount_contacts = amount_contacts;
    }

    public int getActivated_contacts() {
        return activated_contacts;
    }

    public void setActivated_contacts(int activated_contacts) {
        this.activated_contacts = activated_contacts;
    }

    public String getType_name() {
        return type_name;
    }

    public void setType_name(String type_name) {
        this.type_name = type_name;
    }

    public String getActivation_text() {
        return activation_text;
    }

    public void setActivation_text(String activation_text) {
        this.activation_text = activation_text;
    }

    public String getUnsubscribe_url() {
        return unsubscribe_url;
    }

    public void setUnsubscribe_url(String unsubscribe_url) {
        this.unsubscribe_url = unsubscribe_url;
    }

    public ArrayList<ShortUserEmail> getSender() {
        return sender;
    }

    public void setSender(ArrayList<ShortUserEmail> sender) {
        this.sender = sender;
    }

    public int getSender_id() {
        return sender_id;
    }

    public void setSender_id(int sender_id) {
        this.sender_id = sender_id;
    }
    
}

