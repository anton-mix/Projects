package SQLClass.Product;

import MyInterface.ModelInterface;

public class ViewProduct implements ModelInterface{
    int id = -1; 
    int user_id = -1; 
    String product_type = ""; 
    int payment_type = -1;
    String product_name = "";
    String identifier = "";
    Double value = -1.0; 
    
    public ViewProduct(int id, int user_id, String product_type, int payment_type, String product_name, String identifier,double value) {
        this.id = id;
        this.user_id = user_id;
        this.product_type = product_type;
        this.payment_type = payment_type;
        this.product_name = product_name;
        this.identifier = identifier;
        this.value = value;
    }
       
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), product_type, String.valueOf(payment_type), product_name, identifier, String.valueOf(value)};
    }
    
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id),
            String.valueOf(user_id), 
            "'" + product_type.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            String.valueOf(payment_type), 
            "'" + product_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + identifier.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            String.valueOf(value)
        };
    }
     
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "product_type", "payment_type", "product_name", "identifier", "value"}; 
    }

    public String getIdentifier() {
        return identifier;
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    
    
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

    public String getProduct_type() {
        return product_type;
    }

    public void setProduct_type(String product_type) {
        this.product_type = product_type;
    }

    public int getPayment_type() {
        return payment_type;
    }

    public void setPayment_type(int payment_type) {
        this.payment_type = payment_type;
    }

    public String getProduct_name() {
        return product_name;
    }

    public void setProduct_name(String product_name) {
        this.product_name = product_name;
    }

    public Double getValue() {
        return value;
    }

    public void setValue(Double value) {
        this.value = value;
    }
    
    
}
