package SQLClass.Contact;

import MyInterface.ModelInterface;

public class ContactPhone implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int contact_id = -1;
    String phone = "";
    
    public ContactPhone(int id, int user_id, int contact_id, String phone) {
        this.id = id;
        this.user_id = user_id;
        this.contact_id = contact_id;
        this.phone = phone;
    }
    
    public ContactPhone(int user_id, int contact_id, String phone) {
        this.user_id = user_id;
        this.contact_id = contact_id;
        this.phone = phone;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(contact_id), phone};
    }
    
    @Override
    public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id),
            String.valueOf(user_id), 
            String.valueOf(contact_id),
            "'" + phone.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
            };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "contact_id", "phone"}; 
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

    public String getPhone() {
        return phone;
    }

    public void setPhone(String phone) {
        this.phone = phone;
    }
}
