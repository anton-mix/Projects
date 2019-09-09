package SQLClass.Contact;

import MyInterface.ModelInterface;
import java.util.ArrayList;

public class TypeGroupContacts implements ModelInterface {
    int id = -1;
    int user_id = -1; 
    String name = ""; 
    int sender_id = -1;
    int default_group_id = -1;
    String telegram_bot_token = ""; 
    String vk_group_token = "";
    String vk_group_name = "";
    int vk_group_id = -1;
    ArrayList<ShortGroupContacts> groups_contacts;
    ArrayList<ShortUserEmail> sender;
    
    public TypeGroupContacts (int id, int user_id, String name, int sender_id, ArrayList<ShortGroupContacts> groups_contacts, int default_group_id, String telegram_bot_token, String vk_group_token, String vk_group_name, int vk_group_id)  {
        this.id = id;
        this.user_id = user_id;
        this.name = name;   
        this.sender_id = sender_id;
        this.groups_contacts = groups_contacts;
        this.sender = new ArrayList<ShortUserEmail>();
        this.default_group_id = default_group_id;
        this.telegram_bot_token = telegram_bot_token;
        this.vk_group_token = vk_group_token;
        this.vk_group_name = vk_group_name;
        this.vk_group_id = vk_group_id;
    }
    
    public TypeGroupContacts (int user_id, String name, int sender_id, int default_group_id)  {
        this.user_id = user_id;
        this.name = name;  
        this.sender_id = sender_id;
        this.groups_contacts = new ArrayList<ShortGroupContacts>();
        this.sender = new ArrayList<ShortUserEmail>();
        this.default_group_id = default_group_id;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), name, String.valueOf(sender_id), String.valueOf(default_group_id), telegram_bot_token, vk_group_token, vk_group_name, String.valueOf(vk_group_id)};
    }
    
    @Override
    public String[] getSqlArray() {
            return new String[] {
                String.valueOf(id),
                String.valueOf(user_id),
                "'" + name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
                String.valueOf(sender_id),
                String.valueOf(default_group_id)
        };
    }
    
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "name", "sender_id", "default_group_id", "telegram_bot_token", "vk_group_token", "vk_group_name", "vk_group_id"}; 
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

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public ArrayList<ShortGroupContacts> getGroups_contacts() {
        return groups_contacts;
    }

    public void setGroups_contacts(ArrayList<ShortGroupContacts> groups_contacts) {
        this.groups_contacts = groups_contacts;
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

    public int getDefault_group_id() {
        return default_group_id;
    }

    public void setDefault_group_id(int default_group_id) {
        this.default_group_id = default_group_id;
    }

    public String getTelegram_bot_token() {
        return telegram_bot_token;
    }

    public void setTelegram_bot_token(String telegram_bot_token) {
        this.telegram_bot_token = telegram_bot_token;
    }

    public String getVk_group_token() {
        return vk_group_token;
    }

    public void setVk_group_token(String vk_group_token) {
        this.vk_group_token = vk_group_token;
    }

    public String getVk_group_name() {
        return vk_group_name;
    }

    public void setVk_group_name(String vk_group_name) {
        this.vk_group_name = vk_group_name;
    }

    public int getVk_group_id() {
        return vk_group_id;
    }

    public void setVk_group_id(int vk_group_id) {
        this.vk_group_id = vk_group_id;
    }
    
}
