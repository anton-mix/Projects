package SQLClass.Contact;

import MyInterface.ModelInterface;

public class ViewGroupContacts  implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int type_id = -1;
    String type_name = "";
    String group_name = ""; 
    int amount_contacts = -1;
    int activated_contacts = -1;
    
    public ViewGroupContacts(int id, int user_id, int type_id, String group_name, String type_name, int amount_contacts, int activated_contacts){
        this.id = id;
        this.user_id = user_id;
        this.type_id = type_id;
        this.type_name = type_name;
        this.group_name = group_name;
        this.amount_contacts = amount_contacts;
        this.activated_contacts = activated_contacts;
    }
    
    public ViewGroupContacts(int user_id, int type_id, String group_name, String type_name, int amount_contacts, int activated_contacts){
        this.user_id = user_id;
        this.type_id = type_id;
        this.type_name = type_name;
        this.group_name = group_name;
        this.amount_contacts = amount_contacts;
        this.activated_contacts = activated_contacts;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(type_id), group_name, type_name,
            String.valueOf(amount_contacts), String.valueOf(activated_contacts)};
    }
    
    @Override
    public String[] getSqlArray() {
            return new String[] {
                String.valueOf(id),
                String.valueOf(user_id), 
                String.valueOf(type_id), 
                "'" + group_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
                String.valueOf(amount_contacts), 
                String.valueOf(activated_contacts)
        };
    }
    
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "type_id", "group_name", "type_name", "amount_contacts", "activated_contacts"}; 
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
}
