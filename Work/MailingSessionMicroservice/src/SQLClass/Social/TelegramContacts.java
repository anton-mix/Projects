package SQLClass.Social;

import MyInterface.ModelInterface;

public class TelegramContacts implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int contact_id  = -1;
    String first_name = "";
    String last_name = "";
    String username = "";
    int is_subscribed = -1;
    
    public TelegramContacts(int user_id, int contact_id, String first_name, String last_name, String username, int is_subscribed) {
        this.user_id = user_id;
        this.contact_id = contact_id;
        this.first_name = first_name;
        this.last_name = last_name;
        this.username = username;
        this.is_subscribed = is_subscribed;
    } 
    
    public TelegramContacts(int id, int user_id, int contact_id, String first_name, String last_name, String username, int is_subscribed) {
        this.id = id;
        this.user_id = user_id;
        this.contact_id = contact_id;
        this.first_name = first_name;
        this.last_name = last_name;
        this.username = username;
        this.is_subscribed = is_subscribed;
    }  

    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(contact_id), first_name, last_name, username,String.valueOf(is_subscribed)};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id), 
            String.valueOf(contact_id), 
            "'" + first_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + last_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + username.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            String.valueOf(is_subscribed)
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "contact_id", "first_name", "last_name", "username", "username"}; 
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

    public int getContact_id() {
        return contact_id;
    }

    public void setContact_id(int contact_id) {
        this.contact_id = contact_id;
    }

    public String getFirst_name() {
        return first_name;
    }

    public void setFirst_name(String first_name) {
        this.first_name = first_name;
    }

    public String getLast_name() {
        return last_name;
    }

    public void setLast_name(String last_name) {
        this.last_name = last_name;
    }

    public String getUsername() {
        return username;
    }

    public void setUsername(String username) {
        this.username = username;
    }

    public int getIs_subscribed() {
        return is_subscribed;
    }

    public void setIs_subscribed(int is_subscribed) {
        this.is_subscribed = is_subscribed;
    }
}
