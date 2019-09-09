package SQLClass.Contact;

import MyInterface.ModelInterface;

public class ContactSocial implements ModelInterface{
    int id = -1;
    int user_id = -1;
    int contact_id = -1;
    String social = "";
    
    public ContactSocial(int id, int user_id, int contact_id, String social) {
        this.id = id;
        this.user_id = user_id;
        this.contact_id = contact_id;
        this.social = social;
    }
    
    public ContactSocial(int user_id, int contact_id, String social) {
        this.id = id;
        this.user_id = user_id;
        this.contact_id = contact_id;
        this.social = social;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(contact_id), social};
    }
    
    @Override
    public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id),
            String.valueOf(user_id), 
            String.valueOf(contact_id),
            "'" + social.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
            };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "contact_id", "social"}; 
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

    public String getSocial() {
        return social;
    }

    public void setSocial(String social) {
        this.social = social;
    }
}
