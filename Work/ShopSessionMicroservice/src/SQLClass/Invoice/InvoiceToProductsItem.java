package SQLClass.Invoice;

import MyInterface.ModelInterface;
import java.util.List;

public class InvoiceToProductsItem implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int invoice_id = -1;
    int product_id = -1;
    
    
    public InvoiceToProductsItem (int id, int user_id, int invoice_id, int product_id)  {
        this.id = id;
        this.user_id = user_id;
        this.invoice_id = invoice_id;
        this.product_id = product_id;
    }
    
    public InvoiceToProductsItem (int user_id, int invoice_id, int product_id)  {
        this.user_id = user_id;
        this.invoice_id = invoice_id;
        this.product_id = product_id;
    }
    
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(invoice_id), String.valueOf(product_id)};
    }
    
    public String[] getSqlArray() {
        return new String[] {
           String.valueOf(id),
           String.valueOf(user_id),
           String.valueOf(invoice_id),
           String.valueOf(product_id)  
        };
    }
    
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "invoice_id", "product_id"};
    }
@Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }
    
    
    
    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public int getInvoice_id() {
        return invoice_id;
    }

    public void setInvoice_id(int invoice_id) {
        this.invoice_id = invoice_id;
    }

    public int getProduct_id() {
        return product_id;
    }

    public void setProducts_id(int product_id) {
        this.product_id = product_id;
    }
    
    
}
