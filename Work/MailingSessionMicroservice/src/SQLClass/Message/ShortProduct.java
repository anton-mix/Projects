package SQLClass.Message;

public class ShortProduct {
    int id = -1;
    String product_name = "";
    double value = 0;
    
    public ShortProduct() {
        id = -1;
        value = 0;
        product_name = "";
    }
    
    public ShortProduct(int id, double value, String product_name) {
        this.id = id;
        this.value = value;
        this.product_name = product_name;
    }

    public ShortProduct(int id, String product_name) {
        this.id = id;
        this.product_name = product_name;
    }
    
    public String getProduct_name() {
        return product_name;
    }

    public void setProduct_name(String product_name) {
        this.product_name = product_name;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public double getValue() {
        return value;
    }

    public void setValue(double value) {
        this.value = value;
    }
        
}
