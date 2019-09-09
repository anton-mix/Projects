package SQLClass.Product;

import MyInterface.ModelInterface;
import java.util.ArrayList;
import java.util.List;

public class FilterProduct implements ModelInterface {
    int user_id = -1; 
    List<ShortGroupProducts> categories_id; 
    String show_in_catalogue = "";
    String query = "";
    int limit = 0;
    int offset = 0;

    public FilterProduct(int user_id, String show_in_catalogue, String query) {
        this.user_id = user_id;
        categories_id = new ArrayList<ShortGroupProducts>();
        this.show_in_catalogue = show_in_catalogue;
        this.query = query;
    }
    
    public FilterProduct(int user_id, String show_in_catalogue, String query, List<ShortGroupProducts> categories_id, int limit, int offset) {
        this.user_id = user_id;
        categories_id = new ArrayList<ShortGroupProducts>();
        this.show_in_catalogue = show_in_catalogue;
        this.query = query;
        this.categories_id = categories_id;
        this.limit = limit;
        this.offset = offset;
    }
    
    public String getFullText() {
        return "product_name, identifier";
    }
    
    public String[] getArray() {
        return new String[] {String.valueOf(user_id), show_in_catalogue, query, String.valueOf(limit), String.valueOf(offset)};
    }
    
    public String[] getSqlArray() {
        return new String[] {
            String.valueOf(user_id), 
            "'" + show_in_catalogue.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
            };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"user_id", "show_in_catalogue", "query", "limit", "offset"}; 
    }
    
    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public List<ShortGroupProducts> getCategories_id() {
        return categories_id;
    }

    public void setCategories_id(List<ShortGroupProducts> categories_id) {
        this.categories_id = categories_id;
    }

    public String getShow_in_catalogue() {
        return show_in_catalogue;
    }

    public void Show_in_catalogue(String show_in_catalogue) {
        this.show_in_catalogue = show_in_catalogue;
    }

    public String getQuery() {
        return query;
    }

    public void setQuery(String query) {
        this.query = query;
    }

    public int getLimit() {
        return limit;
    }

    public void setLimit(int limit) {
        this.limit = limit;
    }

    public int getOffset() {
        return offset;
    }

    public void setOffset(int offset) {
        this.offset = offset;
    }
    
}
