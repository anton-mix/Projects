package SQLClass.Product;

import MyInterface.ModelInterface;
import java.util.ArrayList;

public class GroupProducts implements ModelInterface {
    int id = -1;
    int user_id = -1;
    String group_name = "";
    ArrayList<ShortProduct> products; 
    
    public GroupProducts(int id, int user_id, String group_name, ArrayList<ShortProduct> products) {
        this.id = id;
        this.user_id = user_id;
        this.group_name = group_name;
        this.products = products;
    }
    
    public GroupProducts(int id, int user_id, String group_name) {
        this.id = id;
        this.user_id = user_id;
        this.group_name = group_name;
        this.products = new ArrayList<ShortProduct>();
    }
    
    public GroupProducts(int user_id, String group_name) {
        this.user_id = user_id;
        this.group_name = group_name;
        this.products = new ArrayList<ShortProduct>();
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), group_name};
    }

    @Override
    public String[] getSqlArray() {
         return new String[] {
            String.valueOf(id),
            String.valueOf(user_id), 
            "'" + group_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            };
    }

    @Override
    public String[] getKeyArray() {
        return new String[] {"id","user_id","group_name"};
    }

    public ArrayList<ShortProduct> getProducts() {
        return products;
    }

    public void setProducts(ArrayList<ShortProduct> products) {
        this.products = products;
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

    public String getGroup_name() {
        return group_name;
    }

    public void setGroup_name(String group_name) {
        this.group_name = group_name;
    }    
}
