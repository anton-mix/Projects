package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.GlobalVariables;
import SQLClass.Product.GroupProducts;
import SQLClass.Product.GroupProductsConverter;
import SQLClass.Product.ShortProduct;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import SQLClass.User.UserSystemInfo;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Connection;
import java.sql.SQLException;

public class GroupsProductsManager implements ManagerInterface {

    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;
    String name_class = "GroupsProductsManager";

    public GroupsProductsManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.global_var = global_var;
        this.connection_with_db = connection_with_db;
        this.log_manager = log_manager;

        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(GroupProducts.class, new GroupProductsConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        gson = builder.create();
    }

    public ResData getGroupsProducts(GroupProducts item) throws MyNotAuthorizedException {
        try {
            GroupProducts[] buff = gson.fromJson(Select("SELECT * FROM GroupsProducts " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson), GroupProducts[].class);
            for (int j = 0; j < buff.length; j++) {
                ShortProduct[] buff2 = gson.fromJson(Select("SELECT id, value, product_name FROM Products WHERE user_id = " + item.getUser_id() + " AND group_id = " + buff[j].getId() + ";", connection_with_db, gson), ShortProduct[].class);
                for (int i = 0; i < buff2.length; i++) {
                    buff[j].getProducts().add(buff2[i]);
                }
            }
            return new ResData(true, gson.toJson(buff));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getGroupsProducts(GroupProducts). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getGroupsProducts(GroupProducts). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getGroupsProducts(GroupProducts). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getGroupsProducts(GroupProducts). Sql select Exception", name_class);
        }
    }

    public ResData getCountGroupsProducts(GroupProducts item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT COUNT(*) as result FROM GroupsProducts " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getCountGroupsProducts(GroupProducts). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getCountGroupsProducts(GroupProducts). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getCountGroupsProducts(GroupProducts). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getCountGroupsProducts(GroupProducts). Sql select Exception", name_class);
        }
    }

    public ResData AddGroupsProducts(GroupProducts item) throws MyNotAuthorizedException {
        String result = "", Buffresult = "";
        try {
            result = Add(item, "GroupsProducts", connection_with_db, gson, true);
            int id = gson.fromJson(result, SQLResult.class).getId();
            for (int i = 0; i < item.getProducts().size(); i++) {
                Buffresult = ExecuteQuery("UPDATE Products SET group_id = " + id + " WHERE user_id = " + item.getUser_id() + " AND id = " + item.getProducts().get(i).getId() + ";", connection_with_db, gson);
            }
            return new ResData(true, result);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddGroupsProducts(GroupProducts). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddGroupsProducts(GroupProducts). Exception in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddGroupsProducts(GroupProducts). Exception", name_class);
        }
    }

    public ResData DeleteGroupsProducts(GroupProducts item) throws MyNotAuthorizedException {
        try {
            UserSystemInfo buff = gson.fromJson(getUserSystemInfo(gson.fromJson("{\"user_id\":\"" + item.getUser_id() + "\"}", UserSystemInfo.class)).getRes(), UserSystemInfo[].class)[0];
            if (buff.getGroup_products_id() == item.getId()) 
                return new ResData(true, gson.toJson(new SQLResult("Ok")));
            
            ExecuteQuery("UPDATE Products SET group_id = " + buff.getGroup_products_id() + " WHERE user_id = " + item.getUser_id() + " AND group_id = " + item.getId() + ";", connection_with_db, gson);
            return new ResData(true, Delete(item, "GroupsProducts", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteGroupsProducts(GroupProducts). Sql delete", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "DeleteGroupsProducts(GroupProducts). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteGroupsProducts(GroupProducts). Sql delete", name_class);
        }
    }

    public ResData UpdateGroupsProducts(GroupProducts item) throws MyNotAuthorizedException {
        try {
            String res = "";
            res = Update(item, "GroupsProducts", connection_with_db, gson);
            for (int i = 0; i < item.getProducts().size(); i++) {
                ExecuteQuery("UPDATE Products SET group_id = " + item.getId() + " WHERE user_id = " + item.getUser_id() + " AND id = " + item.getProducts().get(i).getId() + ";", connection_with_db, gson);
            }
            return new ResData(true, res);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateGroupsProducts(GroupProducts). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateGroupsProducts(GroupProducts). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateGroupsProducts(GroupProducts). Exception", name_class);
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
