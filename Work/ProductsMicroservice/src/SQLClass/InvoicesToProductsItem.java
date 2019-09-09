package SQLClass;

import java.util.List;

public class InvoicesToProductsItem {
    transient int id = -1;
    transient int invoice_id = -1;
    transient int product_id = -1;
    
    enum Invoices_to_Products_enum {id, invoice_id, product_id}
    
    public InvoicesToProductsItem (int id, int invoice_id, int product_id)  {
        this.id = id;
        this.invoice_id = invoice_id;
        this.product_id = product_id;
    }
    
    public InvoicesToProductsItem (int invoice_id, int product_id)  {
        this.invoice_id = invoice_id;
        this.product_id = product_id;
    }
    
    public InvoicesToProductsItem (String json) {
        List<String[]> jsonString = GlobalVariables.Parse_json(json);

        for(int i = 0; i < jsonString.size(); i ++) { 
            
            Invoices_to_Products_enum invoices_to_Products_enum = Invoices_to_Products_enum.valueOf(jsonString.get(i)[0]);
            
            switch(invoices_to_Products_enum) {
                case id:
                    this.id = Integer.valueOf(jsonString.get(i)[1]);
                    break;
                case invoice_id:
                    this.invoice_id = Integer.valueOf(jsonString.get(i)[1]);
                    break;
                case product_id:
                    this.product_id = Integer.valueOf(jsonString.get(i)[1]);
                    break;
            }
        }       
    }
    
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(invoice_id), String.valueOf(product_id)};
    }
    
    public String[] getSqlArray() {
        return new String[] {String.valueOf(id),
           String.valueOf(invoice_id),
           String.valueOf(product_id)  };
    }
    
    public String getSqlStringArray() {
        return "null ," +
            String.valueOf(invoice_id) + "," +
            String.valueOf(product_id) + ",";
    }
    
    public String[] getKeyArray() {
        return new String[] {"id", "invoice_id", "product_id"};
    }
    
    public String getKeyString() {
        return "id, invoice_id, product_id"; 
    }
    
    public String getJson() {
        String[] buffItemString = getArray();
        String[] buffItemKeyString = getKeyArray();
        String Result = "{";
        
        for(int i = 0;i < buffItemString.length;i++) {
            Result = Result +  "\"" + buffItemKeyString[i] + "\":\"" + buffItemString[i] + "\"";
            if(i != buffItemString.length-1)
                    Result = Result + ",";
            }
        return Result +  "}";
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
