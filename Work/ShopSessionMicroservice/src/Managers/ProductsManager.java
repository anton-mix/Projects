package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.GlobalVariables;
import SQLClass.Product.FilterProduct;
import SQLClass.Product.FilterProductConverter;
import SQLClass.Product.ProductModel;
import SQLClass.Product.ProductModelConverter;
import SQLClass.Product.ViewProduct;
import SQLClass.Product.ViewProductConverter;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import SQLClass.User.UserSystemInfo;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;

import java.sql.Connection;
import java.sql.SQLException;

public class ProductsManager implements ManagerInterface {

    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;
    
    String name_class = "ProductsManager";
    
    public ProductsManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.global_var = global_var;
        this.connection_with_db = connection_with_db;
        this.log_manager = log_manager;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(ProductModel.class, new ProductModelConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        builder.registerTypeAdapter(ViewProduct.class, new ViewProductConverter());
        builder.registerTypeAdapter(FilterProduct.class, new FilterProductConverter());
        gson = builder.create();
    }

    public ResData getProducts() {
        try {
            return new ResData(true, Select("SELECT * FROM Products;", connection_with_db, gson));   
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getProducts(). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getProducts(). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getProducts(). Exception", name_class);
        }
    }

    public ResData getProduct(ViewProduct item) throws MyNotAuthorizedException {
         try {
            return new ResData(true, Select("SELECT * FROM Products " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getProducts(ProductFilter). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getProducts(ProductFilter). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getProducts(ProductFilter). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getProducts(ProductFilter). Exception", name_class);
        }
    }
    
    public ResData getViewProducts(FilterProduct item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT * FROM Products " + SqlFilter(item, item.getLimit(), item.getOffset()) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getViewProducts(FilterProduct). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getViewProducts(FilterProduct). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getViewProducts(FilterProduct). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getViewProducts(FilterProduct). Exception", name_class);
        }
    }
    
    public ResData getShortProducts(FilterProduct item) throws MyNotAuthorizedException {
        try {
            ViewProduct[] buf = gson.fromJson(Select("SELECT * FROM Products " + SqlFilter(item, item.getLimit(), item.getOffset()) + ";", connection_with_db, gson) , ViewProduct[].class);
            
            JsonArray products = new JsonArray();
            for(int i = 0; i < buf.length;i++) {
                JsonObject product = new JsonObject();
                product.addProperty("id", buf[i].getId());
                product.addProperty("product_name", buf[i].getProduct_name());
                products.add(product);
            }
            return new ResData(true, gson.toJson(products));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getShortProducts(FilterProduct). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getShortProducts(FilterProduct). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getShortProducts(FilterProduct). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getShortProducts(FilterProduct). Exception", name_class);
        }
    }
    
    public ResData getCountProducts(FilterProduct item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT COUNT(*) as result FROM Products " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getCountProducts(ProductFilter). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getCountProducts(ProductFilter). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getCountProducts(ProductFilter). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getCountProducts(ProductFilter). Exception", name_class);
        }
    }
    
    //Add
    
    public ResData AddProduct(ProductModel item) {
        try {
            UserSystemInfo buff = gson.fromJson(getUserSystemInfo(new UserSystemInfo(item.getUser_id())).getRes(), UserSystemInfo[].class)[0];
            item.setReminder_payment_id(buff.getReminder_payment_id_1());
            item.setReminder_surcharge_id(buff.getReminder_surcharge_id_1());
            return new ResData(true, Add(item, "Products", connection_with_db, gson, true));
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddProduct(ProductModel). Exeption in add products", name_class);
        }
    }
    
    //Delete
    
    public ResData DeleteProduct(ProductModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Delete(item, "Products", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteProduct(ProductModel). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteProduct(ProductModel). Exeption", name_class);
        }
    }
    
    //Update
    
    public ResData UpdateProduct(ProductModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Update(item, "Products", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateProduct(ProductModel). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateProduct(ProductModel). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateProduct(ProductModel). Exeption", name_class);
        }
    }
    
    public ResData Product_select_with_count(FilterProduct item){
        try{
            ResData res, answer;
            
            answer = getViewProducts(item);
            if(!answer.getIs_success()) {
                return answer;
            }
            res = getCountProducts(item);
            if(!res.getIs_success()) {
                return res;
            }
            
            answer.setRes("[" + answer.getRes() + "," + res.getRes().substring(1, res.getRes().length()-1)  + "]");
            answer.setIs_success(res.getIs_success());
            
            return answer;
        }
        catch(MyNotAuthorizedException ex){
            return log_manager.CallError(ex.toString(), "Product_select_with_count(FilterProduct). MyNotAuthorizedException", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "Product_select_with_count(FilterProduct). Exeption", name_class);
        }
    }

    public ResData getUserSystemInfo(UserSystemInfo item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT * FROM UserSystemInfo " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getUserSystemInfo(UserSystemInfo). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getUserSystemInfo(UserSystemInfo). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getUserSystemInfo(UserSystemInfo). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getUserSystemInfo(UserSystemInfo). Sql select Exception", name_class);
        }
    }
    
    public Connection getConnection_with_db() {
        return connection_with_db;
    }

    public void setConnection_with_db(Connection connection_with_db) {
        this.connection_with_db = connection_with_db;
    }
    
    public void setErrorHttp(String ErrorHttp) {
        this.log_manager.setErrorHttp(ErrorHttp);
    }
}
