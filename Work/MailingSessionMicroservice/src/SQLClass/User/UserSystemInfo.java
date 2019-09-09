package SQLClass.User;

import MyInterface.ModelInterface;

public class UserSystemInfo implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int group_products_id = -1;
    int type_group_contacts_id = -1;
    
    public UserSystemInfo(int id, int user_id, int group_products_id, int type_group_contacts_id) {
        this.id = id;
        this.user_id = user_id;
        this.group_products_id = group_products_id;
        this.type_group_contacts_id = type_group_contacts_id;
    }
    
    public UserSystemInfo(int user_id, int group_products_id, int type_group_contacts_id) {
        this.user_id = user_id;
        this.group_products_id = group_products_id;
        this.type_group_contacts_id = type_group_contacts_id;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(group_products_id), String.valueOf(type_group_contacts_id)};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id), 
            String.valueOf(group_products_id), 
            String.valueOf(type_group_contacts_id)
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "group_products_id", "type_group_contacts_id"}; 
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

    public int getGroup_products_id() {
        return group_products_id;
    }

    public void setGroup_products_id(int group_products_id) {
        this.group_products_id = group_products_id;
    }

    public int getType_group_contacts_id() {
        return type_group_contacts_id;
    }

    public void setType_group_contacts_id(int type_group_contacts_id) {
        this.type_group_contacts_id = type_group_contacts_id;
    }
    
}
