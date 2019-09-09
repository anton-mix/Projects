package SQLClass.User;

import MyInterface.ModelInterface;

public class UserSystemInfo implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int group_products_id = -1;
    int type_group_contacts_id = -1;
    int reminder_payment_id_1 = -1;
    int reminder_payment_id_2 = -1;
    int reminder_payment_id_3 = -1;
    int reminder_surcharge_id_1  = -1;
    int reminder_surcharge_id_2  = -1;
    int reminder_surcharge_id_3  = -1;
    
    public UserSystemInfo(int id, int user_id, int group_products_id, int type_group_contacts_id, 
            int reminder_payment_id_1, int reminder_payment_id_2, int reminder_payment_id_3, 
            int reminder_surcharge_id_1, int reminder_surcharge_id_2, int reminder_surcharge_id_3) {
        this.id = id;
        this.user_id = user_id;
        this.group_products_id = group_products_id;
        this.type_group_contacts_id = type_group_contacts_id;
        this.reminder_payment_id_1 = reminder_payment_id_1;
        this.reminder_payment_id_2 = reminder_payment_id_2;
        this.reminder_payment_id_3 = reminder_payment_id_3;
        this.reminder_surcharge_id_1 = reminder_surcharge_id_1;
        this.reminder_surcharge_id_2 = reminder_surcharge_id_2;
        this.reminder_surcharge_id_3 = reminder_surcharge_id_3;
    }
    
    public UserSystemInfo(int user_id, int group_products_id, int type_group_contacts_id,
            int reminder_payment_id_1, int reminder_payment_id_2, int reminder_payment_id_3, 
            int reminder_surcharge_id_1, int reminder_surcharge_id_2, int reminder_surcharge_id_3) {
        this.user_id = user_id;
        this.group_products_id = group_products_id;
        this.type_group_contacts_id = type_group_contacts_id;
        this.reminder_payment_id_1 = reminder_payment_id_1;
        this.reminder_payment_id_2 = reminder_payment_id_2;
        this.reminder_payment_id_3 = reminder_payment_id_3;
        this.reminder_surcharge_id_1 = reminder_surcharge_id_1;
        this.reminder_surcharge_id_2 = reminder_surcharge_id_2;
        this.reminder_surcharge_id_3 = reminder_surcharge_id_3;
    }
    
    public UserSystemInfo(int user_id) {
        this.user_id = user_id;
        this.group_products_id = -1;
        this.type_group_contacts_id = -1;
        this.reminder_payment_id_1 = -1;
        this.reminder_payment_id_2 = -1;
        this.reminder_payment_id_3 = -1;
        this.reminder_surcharge_id_1 = -1;
        this.reminder_surcharge_id_2 = -1;
        this.reminder_surcharge_id_3 = -1;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(group_products_id), String.valueOf(type_group_contacts_id),
        String.valueOf(reminder_payment_id_1), String.valueOf(reminder_payment_id_2), String.valueOf(reminder_payment_id_3),
        String.valueOf(reminder_surcharge_id_1), String.valueOf(reminder_surcharge_id_2), String.valueOf(reminder_surcharge_id_3)};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id), 
            String.valueOf(group_products_id), 
            String.valueOf(type_group_contacts_id),
            String.valueOf(reminder_payment_id_1), 
            String.valueOf(reminder_payment_id_2), 
            String.valueOf(reminder_payment_id_3), 
            String.valueOf(reminder_surcharge_id_1),
            String.valueOf(reminder_surcharge_id_2),
            String.valueOf(reminder_surcharge_id_3)
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "group_products_id", "type_group_contacts_id", 
            "reminder_payment_id_1", "reminder_payment_id_2", "reminder_payment_id_3",
            "reminder_surcharge_id_1", "reminder_surcharge_id_2", "reminder_surcharge_id_3"}; 
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

    public int getReminder_payment_id_1() {
        return reminder_payment_id_1;
    }

    public void setReminder_payment_id_1(int reminder_payment_id_1) {
        this.reminder_payment_id_1 = reminder_payment_id_1;
    }

    public int getReminder_payment_id_2() {
        return reminder_payment_id_2;
    }

    public void setReminder_payment_id_2(int reminder_payment_id_2) {
        this.reminder_payment_id_2 = reminder_payment_id_2;
    }

    public int getReminder_payment_id_3() {
        return reminder_payment_id_3;
    }

    public void setReminder_payment_id_3(int reminder_payment_id_3) {
        this.reminder_payment_id_3 = reminder_payment_id_3;
    }

    public int getReminder_surcharge_id_1() {
        return reminder_surcharge_id_1;
    }

    public void setReminder_surcharge_id_1(int reminder_surcharge_id_1) {
        this.reminder_surcharge_id_1 = reminder_surcharge_id_1;
    }

    public int getReminder_surcharge_id_2() {
        return reminder_surcharge_id_2;
    }

    public void setReminder_surcharge_id_2(int reminder_surcharge_id_2) {
        this.reminder_surcharge_id_2 = reminder_surcharge_id_2;
    }

    public int getReminder_surcharge_id_3() {
        return reminder_surcharge_id_3;
    }

    public void setReminder_surcharge_id_3(int reminder_surcharge_id_3) {
        this.reminder_surcharge_id_3 = reminder_surcharge_id_3;
    }
}
