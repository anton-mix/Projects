package Managers;

import SQLClass.GlobalVariables;
import SQLClass.GroupProducts;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

public class GroupsProductsManager {
    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;

    public GroupsProductsManager(GlobalVariables global_var) throws SQLException, InstantiationException, IllegalAccessException, ClassNotFoundException {
        this.global_var = global_var;
        log_manager = new LogManager(global_var.getHost(), global_var.getQueueNameLog(), "GroupsProductsManager");
        try {
            Class.forName("com.mysql.jdbc.Driver").newInstance();
            String connection_string = "jdbc:mysql://" + global_var.getDbIpSql() + "/" + global_var.getDbName() + "?useSSL=false&user=" + global_var.getDbUserName() + "&password=" + global_var.getDbPassword() + "&useUnicode=true&characterEncoding=UTF-8&serverTimezone=UTC";
            connection_with_db = DriverManager.getConnection(connection_string);
        } catch (SQLException | InstantiationException | IllegalAccessException | ClassNotFoundException ex) {
            log_manager.CallError(ex.toString(), "GroupsProductsManager(). Exception in creating connection");
            throw ex;
        } catch (Exception ex) {
            log_manager.CallError(ex.toString(), "GroupsProductsManager(). Exception in creating connection");
            throw ex;
        }
    } 
    
    public String getGroupsProducts(GroupProducts Item, int Limit, int Offset) {
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
            return log_manager.CallError(ex.toString(), "getGroupsProducts(GroupProducts, int, int). Exeption in parse");
        }

        try {
            String response_string = "";
            statement_with_db = connection_with_db.createStatement();
            result_with_db = statement_with_db.executeQuery("SELECT * FROM GroupsProducts " + where_sql + " ;");

            while (result_with_db.next()) {
                GroupProducts buff_item = new GroupProducts(result_with_db.getInt("id"),
                        result_with_db.getString("user_id"),
                        result_with_db.getString("group_name"));
                if (response_string != "") {
                    response_string = response_string + ",";
                }
                response_string = response_string + buff_item.getJson();
            }
            result_with_db.close();

            statement_with_db.close();
            return "[" + response_string + "]";

        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getGroupsProducts(GroupProducts, int, int). Sql all select");
        }
    }
    
    public String getCountGroupsProducts(GroupProducts Item) {
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
            return log_manager.CallError(ex.toString(), "getCountGroupsProducts(GroupProducts, int, int). Exeption in parse");
        }

        try {
            int count = 0;

            StatementWithDB = connection_with_db.createStatement();

            result_with_db = StatementWithDB.executeQuery("SELECT Count(*) FROM GroupsProducts " + WhereSql + " ;");
            result_with_db.next();
            count = result_with_db.getInt(1);
            
            result_with_db.close();
            StatementWithDB.close();
            
            return "[{\"Result\":\" " + count + "\" }]";
        } 
        catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getCountGroupsProducts(GroupProducts, int, int). Sql all select");
        }
    }

    public String AddGroupProducts(GroupProducts item) {
        Statement StatementWithDB = null;
        
        try {
            StatementWithDB = connection_with_db.createStatement();
            StatementWithDB.executeUpdate("INSERT INTO GroupsProducts ( " + item.getKeyString() + " ) VALUES (" + item.getSqlStringArray() + ");");
                
            StatementWithDB.close();
            connection_with_db.close();
                
            return "[{\"Result\":\"Ok\"}]";
        }
        catch(SQLException ex)  {
            return log_manager.CallError(ex.toString(), "AddGroupProducts(GroupProducts). Add sql");
        }
    }
    
    public String DeleteGroupProducts(int id) {
        Statement statement_with_db = null;
        
        try {
            statement_with_db = connection_with_db.createStatement();
            statement_with_db.executeUpdate("DELETE FROM GroupsProducts WHERE id=" + id + ";");

            statement_with_db.close();
            
            return "[{\"Result\":\"Ok\"}]";
        }
        catch(SQLException ex)  {
            return log_manager.CallError(ex.toString(), "DeleteGroupProducts(int). Sql del");
        }
    }
    
    public String UpdateGroupProducts(GroupProducts item){
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
            return log_manager.CallError(ex.toString(), "UpdateGroupProducts(GroupProducts). Exception In parse");
        }
        
        try {
            statement_with_db = connection_with_db.createStatement();
            statement_with_db.executeUpdate("UPDATE GroupsProducts SET " + set_string + " WHERE " + key_product[0] + "=" + attribute_product[0] + " ;");
            statement_with_db.close();            
            return "[{\"Result\":\"Ok\"}]";
        }
        catch(SQLException ex)  {
            return log_manager.CallError(ex.toString(), "UpdateGroupProducts(GroupProducts). Sql delete");
        }
    }
}
