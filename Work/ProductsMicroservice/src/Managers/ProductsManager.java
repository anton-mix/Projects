package Managers;

import SQLClass.Product;
import SQLClass.GlobalVariables;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

public class ProductsManager {

    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;

    public ProductsManager(GlobalVariables global_var) throws SQLException, InstantiationException, IllegalAccessException, ClassNotFoundException {
        this.global_var = global_var;
        log_manager = new LogManager(global_var.getHost(), global_var.getQueueNameLog(), "ProductManager");
        try {
            Class.forName("com.mysql.jdbc.Driver").newInstance();
            String connection_string = "jdbc:mysql://" + global_var.getDbIpSql() + "/" + global_var.getDbName() + "?useSSL=false&user=" + global_var.getDbUserName() + "&password=" + global_var.getDbPassword() + "&useUnicode=true&characterEncoding=UTF-8&serverTimezone=UTC";
            connection_with_db = DriverManager.getConnection(connection_string);
        } catch (SQLException | InstantiationException | IllegalAccessException | ClassNotFoundException ex) {
            log_manager.CallError(ex.toString(), "ProductManager(). Exception in creating connection");
            throw ex;
        } catch (Exception ex) {
            log_manager.CallError(ex.toString(), "ProductManager(). Exception in creating connection");
            throw ex;
        }
    }

    public String getProducts() {
        ResultSet result_with_db = null;
        String response_string = "";
        try {
            Statement statement_with_db = connection_with_db.createStatement();
            result_with_db = statement_with_db.executeQuery("SELECT * FROM Products;");
            while (result_with_db.next()) {
                Product buff_item = new Product(result_with_db.getInt("id"),
                        result_with_db.getString("user_id"),
                        result_with_db.getString("product_type"),
                        result_with_db.getString("payment_type"),
                        result_with_db.getString("order_page"),
                        result_with_db.getString("product_name"),
                        result_with_db.getDouble("value"),
                        result_with_db.getString("some_payments"),
                        result_with_db.getDouble("prepayment_min"),
                        result_with_db.getString("image"),
                        result_with_db.getString("show_in_catalogue"),
                        result_with_db.getString("product_description"),
                        result_with_db.getDouble("amount_of_expenses"),
                        result_with_db.getDouble("percents_of_expenses"),
                        result_with_db.getString("groups_after_prepayment"),
                        result_with_db.getString("message_topic_after_prepayment"),
                        result_with_db.getString("message_type_after_prepayment"),
                        result_with_db.getString("message_after_prepayment"),
                        result_with_db.getString("groups_after_payment"),
                        result_with_db.getString("thanks_page"),
                        result_with_db.getString("message_topic_after_payment"),
                        result_with_db.getString("message_type_after_payment"),
                        result_with_db.getString("message_after_payment"),
                        result_with_db.getString("getting_page"),
                        result_with_db.getString("pincodes"),
                        result_with_db.getString("tax"),
                        result_with_db.getString("vendor_code"),
                        result_with_db.getString("company"),
                        result_with_db.getString("country"),
                        result_with_db.getString("ordering"),
                        result_with_db.getString("bill_type"),
                        result_with_db.getTimestamp("bill_date"),
                        result_with_db.getDouble("comission"),
                        result_with_db.getDouble("comission_percents"),
                        result_with_db.getString("show_for_partners"),
                        result_with_db.getDouble("employee_reward"),
                        result_with_db.getDouble("employee_reward_percents"),
                        result_with_db.getString("employee_instruction"),
                        result_with_db.getString("product_type_for_add"),
                        result_with_db.getString("personal_template"),
                        result_with_db.getString("api_url_order"),
                        result_with_db.getString("api_url_prepayment"),
                        result_with_db.getString("api_url_payment"),
                        result_with_db.getString("api_url_cancel"),
                        result_with_db.getString("api_url_return"),
                        result_with_db.getString("api_url_alert_return"));
                if (response_string != null) {
                    response_string = response_string + ",";
                }
                response_string = response_string + buff_item.getJson();
            }
            
            statement_with_db.close();
            result_with_db.close();
            return "[" + response_string + "]";

        } 
        catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getProduct(). Exception in getProduct");
        }
    }

    public String getProducts(Product Item, int Limit, int Offset) {
        Statement statement_with_db = null;
        ResultSet result_with_db = null;
        String where_sql = "";
        String[] attribute_products = Item.getSqlArray();
        String[] key_products = Item.getKeyArray();

        try {
            for (int i = 0; i < attribute_products.length; i++) {
                if (!attribute_products[i].equals("''") && !attribute_products[i].equals("-1") && !attribute_products[i].equals("-1.0") && !attribute_products[i].equals("'1970-01-02 00:00:00.0'")) {
                    if (where_sql == "") {
                        where_sql = "WHERE " + key_products[i] + "=" + attribute_products[i];
                    } else {
                        where_sql = where_sql + " AND " + key_products[i] + "=" + attribute_products[i] + " ";
                    }
                }
            }
            
            if (Limit > 0) {
                where_sql = where_sql + " LIMIT " + Limit;
                if (Offset >= 0) {
                    where_sql = where_sql + " OFFSET " + Offset;
                }
            }
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getProduct(Product, int, int). Exeption in parse");
        }

        try {
            String response_string = "";
            statement_with_db = connection_with_db.createStatement();
            result_with_db = statement_with_db.executeQuery("SELECT * FROM Products " + where_sql + " ;");

            while (result_with_db.next()) {
                Product buff_item = new Product(result_with_db.getInt("id"),
                        result_with_db.getString("user_id"),
                        result_with_db.getString("product_type"),
                        result_with_db.getString("payment_type"),
                        result_with_db.getString("order_page"),
                        result_with_db.getString("product_name"),
                        result_with_db.getDouble("value"),
                        result_with_db.getString("some_payments"),
                        result_with_db.getDouble("prepayment_min"),
                        result_with_db.getString("image"),
                        result_with_db.getString("show_in_catalogue"),
                        result_with_db.getString("product_description"),
                        result_with_db.getDouble("amount_of_expenses"),
                        result_with_db.getDouble("percents_of_expenses"),
                        result_with_db.getString("groups_after_prepayment"),
                        result_with_db.getString("message_topic_after_prepayment"),
                        result_with_db.getString("message_type_after_prepayment"),
                        result_with_db.getString("message_after_prepayment"),
                        result_with_db.getString("groups_after_payment"),
                        result_with_db.getString("thanks_page"),
                        result_with_db.getString("message_topic_after_payment"),
                        result_with_db.getString("message_type_after_payment"),
                        result_with_db.getString("message_after_payment"),
                        result_with_db.getString("getting_page"),
                        result_with_db.getString("pincodes"),
                        result_with_db.getString("tax"),
                        result_with_db.getString("vendor_code"),
                        result_with_db.getString("company"),
                        result_with_db.getString("country"),
                        result_with_db.getString("ordering"),
                        result_with_db.getString("bill_type"),
                        result_with_db.getTimestamp("bill_date"),
                        result_with_db.getDouble("comission"),
                        result_with_db.getDouble("comission_percents"),
                        result_with_db.getString("show_for_partners"),
                        result_with_db.getDouble("employee_reward"),
                        result_with_db.getDouble("employee_reward_percents"),
                        result_with_db.getString("employee_instruction"),
                        result_with_db.getString("product_type_for_add"),
                        result_with_db.getString("personal_template"),
                        result_with_db.getString("api_url_order"),
                        result_with_db.getString("api_url_prepayment"),
                        result_with_db.getString("api_url_payment"),
                        result_with_db.getString("api_url_cancel"),
                        result_with_db.getString("api_url_return"),
                        result_with_db.getString("api_url_alert_return"));
                if (response_string != "") {
                    response_string = response_string + ",";
                }
                response_string = response_string + buff_item.getJson();
            }
            result_with_db.close();

            statement_with_db.close();
            return "[" + response_string + "]";

        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getProduct(Product, int, int). Sql all select");
        }
    }
    
    public String searchProducts (String SearchString, Product Item, int Limit, int Offset) {
        Statement statement_with_db = null;
        ResultSet result_with_db = null;
        String where_sql = "";
        String[] attribute_products = Item.getSqlArray();
        String[] key_products = Item.getKeyArray();

        try {
            for (int i = 0; i < attribute_products.length; i++) {
                if (!attribute_products[i].equals("''") && !attribute_products[i].equals("-1") && !attribute_products[i].equals("-1.0") && !attribute_products[i].equals("'1970-01-02 00:00:00.0'")) {
                    if (where_sql == "") {
                        where_sql = "WHERE " + key_products[i] + "=" + attribute_products[i];
                    } else {
                        where_sql = where_sql + " AND " + key_products[i] + "=" + attribute_products[i] + " ";
                    }
                }
            }
            
            if (Limit > 0) {
                where_sql = where_sql + " LIMIT " + Limit;
                if (Offset >= 0) {
                    where_sql = where_sql + " OFFSET " + Offset;
                }
            }
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getProduct(Product, int, int). Exeption in parse");
        }

        try {
            String response_string = "";
            statement_with_db = connection_with_db.createStatement();
            result_with_db = statement_with_db.executeQuery("SELECT * FROM Products " + where_sql + " ;");

            while (result_with_db.next()) {
                Product buff_item = new Product(result_with_db.getInt("id"),
                        result_with_db.getString("user_id"),
                        result_with_db.getString("product_type"),
                        result_with_db.getString("payment_type"),
                        result_with_db.getString("order_page"),
                        result_with_db.getString("product_name"),
                        result_with_db.getDouble("value"),
                        result_with_db.getString("some_payments"),
                        result_with_db.getDouble("prepayment_min"),
                        result_with_db.getString("image"),
                        result_with_db.getString("show_in_catalogue"),
                        result_with_db.getString("product_description"),
                        result_with_db.getDouble("amount_of_expenses"),
                        result_with_db.getDouble("percents_of_expenses"),
                        result_with_db.getString("groups_after_prepayment"),
                        result_with_db.getString("message_topic_after_prepayment"),
                        result_with_db.getString("message_type_after_prepayment"),
                        result_with_db.getString("message_after_prepayment"),
                        result_with_db.getString("groups_after_payment"),
                        result_with_db.getString("thanks_page"),
                        result_with_db.getString("message_topic_after_payment"),
                        result_with_db.getString("message_type_after_payment"),
                        result_with_db.getString("message_after_payment"),
                        result_with_db.getString("getting_page"),
                        result_with_db.getString("pincodes"),
                        result_with_db.getString("tax"),
                        result_with_db.getString("vendor_code"),
                        result_with_db.getString("company"),
                        result_with_db.getString("country"),
                        result_with_db.getString("ordering"),
                        result_with_db.getString("bill_type"),
                        result_with_db.getTimestamp("bill_date"),
                        result_with_db.getDouble("comission"),
                        result_with_db.getDouble("comission_percents"),
                        result_with_db.getString("show_for_partners"),
                        result_with_db.getDouble("employee_reward"),
                        result_with_db.getDouble("employee_reward_percents"),
                        result_with_db.getString("employee_instruction"),
                        result_with_db.getString("product_type_for_add"),
                        result_with_db.getString("personal_template"),
                        result_with_db.getString("api_url_order"),
                        result_with_db.getString("api_url_prepayment"),
                        result_with_db.getString("api_url_payment"),
                        result_with_db.getString("api_url_cancel"),
                        result_with_db.getString("api_url_return"),
                        result_with_db.getString("api_url_alert_return"));
                if (response_string != "") {
                    response_string = response_string + ",";
                }
                response_string = response_string + buff_item.getJson();
            }
            result_with_db.close();

            statement_with_db.close();
            return "[" + response_string + "]";

        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getProduct(Product, int, int). Sql all select");
        }
    }
    
    
    public String getProductsByGroupId(int id) {
        Statement statement_with_db = null;
        ResultSet result_with_db = null;
     
        try {
            String response_string = "";
            statement_with_db = connection_with_db.createStatement();
            result_with_db = statement_with_db.executeQuery("SELECT a.id, a.product_name FROM Products a, Products_to_GroupProducts b, GroupsProducts c WHERE (c.id = " + id + ") and (c.id = b.groupProducts_id) and (a.id = b.product_id);");

            while (result_with_db.next()) {
                Product buff_item = new Product("{ \"id\":" + result_with_db.getInt("id") + ",\"product_name\":\"" + result_with_db.getString("product_name") + "\"}");
                if (response_string != "") {
                    response_string = response_string + ",";
                }
                response_string = response_string + buff_item.getJson();
            }
            result_with_db.close();

            statement_with_db.close();
            return "[" + response_string + "]";

        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getProduct(Product, int, int). Sql all select");
        }
    }
    

    public String getCountProducts(Product Item) {
        Statement StatementWithDB = null;
        ResultSet result_with_db = null;
        String WhereSql = "";
        String[] AttributeProduct = Item.getSqlArray();
        String[] KeyProduct = Item.getKeyArray();

        try {
            for (int i = 0; i < AttributeProduct.length; i++) {
                if (!AttributeProduct[i].equals("''") && !AttributeProduct[i].equals("-1") && !AttributeProduct[i].equals("-1.0") && !AttributeProduct[i].equals("'1970-01-02 00:00:00.0'")) {
                    if (WhereSql == "") {
                        WhereSql = "WHERE " + KeyProduct[i] + "=" + AttributeProduct[i];
                    } else {
                        WhereSql = WhereSql + " AND " + KeyProduct[i] + "=" + AttributeProduct[i] + " ";
                    }
                }
            }
        } 
        catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getCountProduct(Product, int, int). Exeption in parse" + WhereSql);
        }

        try {
            int count = 0;

            StatementWithDB = connection_with_db.createStatement();

            result_with_db = StatementWithDB.executeQuery("SELECT Count(*) FROM Products " + WhereSql + " ;");
            result_with_db.next();
            count = result_with_db.getInt(1);
            
            result_with_db.close();
            StatementWithDB.close();
            
            return "[{\"Result\":\" " + count + "\" }]";
        } 
        catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getProduct(Product, int, int). Sql all select");
        }
    }

    public String AddProduct(Product item) {
        Statement StatementWithDB = null;
        
        String[] AttributeProduct = item.getSqlArray();
        String[] KeyProduct = item.getKeyArray();
        String NeedAttributeProduct = "", NeedKeyProduct = "";
        try {
            for (int i = 0; i < AttributeProduct.length; i++) {
                if (!AttributeProduct[i].equals("''") && !AttributeProduct[i].equals("-1") && !AttributeProduct[i].equals("-1.0") && !AttributeProduct[i].equals("'1970-01-02 00:00:00.0'")) {
                    if (NeedAttributeProduct == "") {
                        NeedAttributeProduct = AttributeProduct[i];
                        NeedKeyProduct = KeyProduct[i];
                    } else {
                        NeedAttributeProduct = NeedAttributeProduct + "," +  AttributeProduct[i];
                        NeedKeyProduct = NeedKeyProduct + "," + KeyProduct[i];
                    }
                }
            }
        } 
        catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddProduct(Product). Exeption in parse");
        }
        
        try {
            StatementWithDB = connection_with_db.createStatement();
            StatementWithDB.executeUpdate("INSERT INTO Products ( " + NeedKeyProduct + " ) VALUES (" + NeedAttributeProduct + ");");
                
            StatementWithDB.close();
                
            return "[{\"Result\":\"Ok\"}]";
        }
        catch(SQLException ex)  {
            return log_manager.CallError(ex.toString(), "AddProduct(Product). Add sql");
        }
    }
    
    public String DeleteProduct(int id) {
        Statement statement_with_db = null;
        
        try {
            statement_with_db = connection_with_db.createStatement();
            statement_with_db.executeUpdate("DELETE FROM Products WHERE id=" + id + ";");

            statement_with_db.close();
            
            return "[{\"Result\":\"Ok\"}]";
        }
        catch(SQLException ex)  {
            return log_manager.CallError(ex.toString(), "DeleteProduct(Product). Sql del");
        }
    }
    
    public String UpdateProduct(Product item){
        Statement statement_with_db = null;
        String set_string = "" ;
        String[] attribute_product = item.getSqlArray();
        String[] key_product = item.getKeyArray();
        
        try {
            for (int i = 1; i < attribute_product.length; i++) {
                if (!attribute_product[i].equals("''") && !attribute_product[i].equals("-1") && !attribute_product[i].equals("-1.0") && !attribute_product[i].equals("'1970-01-02 00:00:00.0'")) {
                    if (set_string == "") {
                        set_string = set_string + key_product[i] + "=" + attribute_product[i];
                    } else {
                        set_string = set_string + ", " + key_product[i] + "=" + attribute_product[i] + " ";
                    }
                }
            }
        } 
        catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateProduct. Exception In parse");
        }
        
        try {
            statement_with_db = connection_with_db.createStatement();
            statement_with_db.executeUpdate("UPDATE Products SET " + set_string + " WHERE " + key_product[0] + "=" + attribute_product[0] + " ;");
            statement_with_db.close();            
            return "[{\"Result\":\"Ok\"}]";
        }
        catch(SQLException ex)  {
            return log_manager.CallError(ex.toString(), "UpdateProduct(Product). Sql delete");
        }
    }
}
