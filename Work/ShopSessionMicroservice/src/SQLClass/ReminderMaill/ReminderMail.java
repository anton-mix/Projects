package SQLClass.ReminderMaill;

import MyInterface.ModelInterface;
import java.sql.Timestamp;

public class ReminderMail implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int contact_id = -1;
    int maill_tamplate_id = -1;
    int product_id = -1;
    Timestamp date_send;
    
    public ReminderMail(int id, int user_id, int contact_id, int maill_tamplate_id, int product_id, Timestamp date_send) {
        this.id = id;
        this.user_id = user_id;
        this.contact_id = contact_id;
        this.maill_tamplate_id = maill_tamplate_id;
        this.product_id = product_id;
        this.date_send = date_send;
    }
    
    public ReminderMail(int user_id, int contact_id, int maill_tamplate_id, int product_id, Timestamp date_send) {
        this.user_id = user_id;
        this.contact_id = contact_id;
        this.maill_tamplate_id = maill_tamplate_id;
        this.product_id = product_id;
        this.date_send = date_send;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(contact_id), String.valueOf(maill_tamplate_id),
            String.valueOf(product_id), date_send.toString()};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id),
            String.valueOf(contact_id), 
            String.valueOf(maill_tamplate_id), 
            String.valueOf(product_id), 
            date_send.toString(),
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "contact_id", "maill_tamplate_id", "product_id", "date_send"}; 
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

    public int getMaill_tamplate_id() {
        return maill_tamplate_id;
    }

    public void setMaill_tamplate_id(int maill_tamplate_id) {
        this.maill_tamplate_id = maill_tamplate_id;
    }

    public int getProduct_id() {
        return product_id;
    }

    public void setProduct_id(int product_id) {
        this.product_id = product_id;
    }

    public Timestamp getDate_send() {
        return date_send;
    }

    public void setDate_send(Timestamp data_send) {
        this.date_send = data_send;
    }
}
