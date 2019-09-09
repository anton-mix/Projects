package SQLClass;

import java.util.List;

public class ProductsToGroupProductsItem {
    transient int id = -1;
    transient int product_id = -1;
    transient int groupProducts_id = -1;
    
    enum Products_to_GroupProducts_enum {id, product_id, groupProducts_id}
    
    public ProductsToGroupProductsItem (int id, int product_id, int groupProducts_id)  {
        this.id = id;
        this.product_id = product_id;
        this.groupProducts_id = groupProducts_id;
    }
    
    public ProductsToGroupProductsItem (int product_id, int groupProducts_id)  {
        this.product_id = product_id;
        this.groupProducts_id = groupProducts_id;
    }
    
    public ProductsToGroupProductsItem (String json) {
        
        List<String[]> jsonString = GlobalVariables.Parse_json(json);

        for(int i = 0; i < jsonString.size(); i ++) {

            Products_to_GroupProducts_enum products_to_GroupProducts_enum = Products_to_GroupProducts_enum.valueOf(jsonString.get(i)[0]);
            
            switch(products_to_GroupProducts_enum) {
                case id:
                    this.id = Integer.valueOf(jsonString.get(i)[1]);
                    break;
                case product_id:
                    this.product_id = Integer.valueOf(jsonString.get(i)[1]);
                    break;
                case groupProducts_id:
                    this.groupProducts_id = Integer.valueOf(jsonString.get(i)[1]);
                    break;
            }
        }       
    }
    
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(product_id), String.valueOf(groupProducts_id)};
    }
    
    public String[] getSqlArray() {
        return new String[] {String.valueOf(id),
           String.valueOf(product_id),
           String.valueOf(groupProducts_id)  };
    }
    
    public String getSqlStringArray() {
        return "null ," +
            String.valueOf(product_id) + "," +
            String.valueOf(groupProducts_id) + ",";
    }
    
    public String[] getKeyArray() {
        return new String[] {"id", "products_id", "groupProducts_id"};
    }
    
    public String getKeyString() {
        return "id, products_id, groupProducts_id"; 
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

    public int getProduct_id() {
        return product_id;
    }

    public void setProduct_id(int product_id) {
        this.product_id = product_id;
    }

    public int getGroupProducts_id() {
        return groupProducts_id;
    }

    public void setGroupProducts_id(int groupProducts_id) {
        this.groupProducts_id = groupProducts_id;
    }
}
