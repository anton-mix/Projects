package SQLClass.ReminderSurcharge;

import MyInterface.ModelInterface;
import SQLClass.Product.ShortProduct;
import SQLClass.ReminderMaill.ReminderMessageView;
import java.util.ArrayList;

public class ReminderSurchargeModel implements ModelInterface {
    int id = -1;
    int user_id = -1;
    String name = "";
    int count;
    int limit = -1;
    int offset = -1;
    
    ArrayList<ShortProduct> products;
    ArrayList<ReminderMessageView> messages;
    
    public ReminderSurchargeModel(int id, int user_id, String name, int count, ArrayList<ShortProduct> products, ArrayList<ReminderMessageView> messages) {
        this.id = id;
        this.user_id = user_id;
        this.name = name;
        this.count = count;
        this.products = products;
        this.messages = messages;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), name, String.valueOf(count)};
    }

    @Override
    public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id),
            String.valueOf(user_id),
            "'" + name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            String.valueOf(count),
        };
    }
        
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "name", "count"};
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

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public int getCount() {
        return count;
    }

    public void setCount(int count) {
        this.count = count;
    }

    public ArrayList<ShortProduct> getProducts() {
        return products;
    }

    public void setProducts(ArrayList<ShortProduct> products) {
        this.products = products;
    }

    @Override
    public int getLimit() {
        return limit;
    }

    public void setLimit(int limit) {
        this.limit = limit;
    }

    @Override
    public int getOffset() {
        return offset;
    }

    public void setOffset(int offset) {
        this.offset = offset;
    }

    public ArrayList<ReminderMessageView> getMessages() {
        return messages;
    }

    public void setMessages(ArrayList<ReminderMessageView> messages) {
        this.messages = messages;
    }
}
