
package SQLClass.Contact;

import MyInterface.ModelInterface;

public class ContactsToGroupContact implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int contacts_id = -1; 
    int group_contacts_id = -1;
    
    public ContactsToGroupContact (int id, int user_id, int contacts_id, int group_contacts_id){
        this.id = id;
        this.user_id = user_id;
        this.contacts_id = contacts_id;
        this.group_contacts_id = group_contacts_id;
    }
    
    public ContactsToGroupContact (int user_id, int contacts_id, int group_contacts_id) {
        this.user_id = user_id;
        this.contacts_id = contacts_id;
        this.group_contacts_id = group_contacts_id;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(contacts_id), String.valueOf(group_contacts_id)};
    }
    
    @Override
    public String[] getSqlArray() {
            return new String[] {
                String.valueOf(id), 
                String.valueOf(user_id),
                String.valueOf(contacts_id), 
                String.valueOf(group_contacts_id)         
        };
    }
    
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "contacts_id","group_contacts_id"}; 
    }

    @Override
    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public int getContacts_id() {
        return contacts_id;
    }

    public void setContacts_id(int contacts_id) {
        this.contacts_id = contacts_id;
    }

    public int getGroup_contacts_id() {
        return group_contacts_id;
    }

    public void setGroup_contacts_id(int group_contacts_id) {
        this.group_contacts_id = group_contacts_id;
    }

    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }
}
