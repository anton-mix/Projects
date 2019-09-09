package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.GlobalVariables;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import SQLClass.User.UserSystemInfo;
import SQLClass.User.UserSystemInfoConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Connection;
import java.sql.SQLException;

public class UserSystemInfoManager implements ManagerInterface {
    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;
    
    String name_class = "UserSystemInfoManager";
    
    public UserSystemInfoManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.global_var = global_var;
        this.connection_with_db = connection_with_db;
        this.log_manager = log_manager;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(UserSystemInfo.class, new UserSystemInfoConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        gson = builder.create();
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
            return log_manager.CallError(ex.toString(), "getUserSystemInfo(UserSystemInfo). Exception", name_class);
        }
    }
    
    //Add
    
    public ResData addUserSystemInfo(UserSystemInfo item) {
        try {
            return new ResData(false, Add(item, "UserSystemInfo", connection_with_db, gson, true));
        }  catch (MyNotAuthorizedException | MyParseException | SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddUserSystemInfo(UserSystemInfo). Exeption in add products", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddUserSystemInfo(UserSystemInfo). Exception", name_class);
        }
    }
    
    //Delete
    
    public ResData DeleteUserSystemInfo(UserSystemInfo item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Delete(item, "UserSystemInfo", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteUserSystemInfo(UserSystemInfo). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteUserSystemInfo(UserSystemInfo). Exception", name_class);
        }
    }
    
    //Update
    
    public ResData UpdateUserSystemInfo(UserSystemInfo item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Update(item, "UserSystemInfo", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateProduct(ProductModel). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateProduct(ProductModel). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateProduct(ProductModel). Exception", name_class);
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
