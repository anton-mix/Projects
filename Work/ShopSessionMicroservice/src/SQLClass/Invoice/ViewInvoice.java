package SQLClass.Invoice;

import MyInterface.ModelInterface;
import SQLClass.Product.ShortProduct;
import java.sql.Timestamp;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

public class ViewInvoice implements ModelInterface {
    int id = -1;
    String email = "";
    String name = ""; 
    String phone = "";
    int user_id = -1;
    double amount = -1.0; 
    Timestamp date_create; 
    String payment_status = "";
    ArrayList<ShortProduct> products = new ArrayList<ShortProduct>();
      
    DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    Date date = dateFormat.parse("1970-01-02 00:00:00");
    long time = date.getTime();
    
    public ViewInvoice(int id, int user_id, String email, String name, String phone, double amount, Timestamp date_create, 
             String payment_status) throws ParseException {
        this.id = id;
        this.user_id = user_id;
        this.email = email;
        this.name = name;
        this.phone = phone;
        this.amount = amount;
        if(date_create == null) 
            this.date_create = new Timestamp(time);
        else
            this.date_create = date_create;
        this.payment_status = payment_status;
        products = new ArrayList<ShortProduct>();
    }
    
    public ViewInvoice(int id, int user_id, String email, String name, String phone, double amount, Timestamp date_create, 
            String payment_status, ArrayList<ShortProduct> products) throws ParseException {
        this.id = id;
        this.user_id = user_id;
        this.email = email;
        this.name = name;
        this.phone = phone;
        this.amount = amount;
        if(date_create == null) 
            this.date_create = new Timestamp(time);
        else
            this.date_create = date_create;
        this.payment_status = payment_status;
        this.products = products;
    }
       
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), email, name, phone, String.valueOf(amount), date_create.toString(), payment_status};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id),
            String.valueOf(user_id), 
            "'" + email.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + phone.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            String.valueOf(amount), 
            "'" + date_create.toString() + "'", 
            "'" + payment_status.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "email", "name", "phone", "amount", "date_create", "payment_status"}; 
    }
    
    @Override
    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public String getEmail() {
        return email;
    }

    public void setEmail(String email) {
        this.email = email;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getPhone() {
        return phone;
    }

    public void setPhone(String phone) {
        this.phone = phone;
    }

    public String getPayment_status() {
        return payment_status;
    }

    public void setPayment_status(String payment_status) {
        this.payment_status = payment_status;
    }

    public Timestamp getDate_create() {
        return date_create;
    }

    public void setDate_create(Timestamp date_create) {
        this.date_create = date_create;
    }

    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public ArrayList<ShortProduct> getProducts() {
        return products;
    }

    public void setProducts(ArrayList<ShortProduct> products) {
        this.products = products;
    }

    public double getAmount() {
        return amount;
    }

    public void setAmount(double amount) {
        this.amount = amount;
    }
}
