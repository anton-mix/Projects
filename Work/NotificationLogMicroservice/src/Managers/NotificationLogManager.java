package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.GlobalVariables;
import SQLClass.NotificationLog.NotificationLogModel;
import SQLClass.NotificationLog.NotificationLogModelConverter;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Connection;
import java.sql.SQLException;

public class NotificationLogManager implements ManagerInterface {
    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;
    String name_class = "NotificationLogManager";

    public NotificationLogManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.log_manager = log_manager;
        this.connection_with_db = connection_with_db;
        this.global_var = global_var;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(NotificationLogModel.class, new NotificationLogModelConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        gson = builder.create();
    }

    public ResData getNotificationLog(NotificationLogModel item) throws MyNotAuthorizedException {
        try {           
            return new ResData(true, Select("SELECT * FROM NotificationLogs;", connection_with_db, gson)); 
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getNotificationLog(NotificationLogModel). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getNotificationLog(NotificationLogModel). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getNotificationLog(NotificationLogModel). Exception", name_class);
        }
    }
    
    public ResData getCountNotificationLog(NotificationLogModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT COUNT(*) as result FROM NotificationLogs " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getCountNotificationLog(NotificationLogModel). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getCountNotificationLog(NotificationLogModel). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getCountNotificationLog(NotificationLogModel). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getCountNotificationLog(NotificationLogModel). Exception", name_class);
        }
    }
    
    //Add
    
    public ResData AddNotificationLog(NotificationLogModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Add(item, "NotificationLogs", connection_with_db, gson, true));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddNotificationLog(NotificationLogModel). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddNotificationLog(NotificationLogModel). Exeption in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddNotificationLog(NotificationLogModel). Exception", name_class);
        }
    }
    
    //Delete
    
    public ResData DeleteNotificationLog(NotificationLogModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Delete(item, "NotificationLogs", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteNotificationLog(NotificationLogModel). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteNotificationLog(NotificationLogModel). Exception", name_class);
        }
    }
    
     //Update
    
    public ResData UpdateNotificationLog(NotificationLogModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Update(item, "NotificationLogs", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateNotificationLog(NotificationLogModel). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateNotificationLog(NotificationLogModel). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateNotificationLog(NotificationLogModel). Exception", name_class);
        }
    }
    
    public ResData Notification_log_select_with_count(NotificationLogModel item) throws MyNotAuthorizedException {
        ResData answer, res;
        
        answer = getNotificationLog(item);
        if(!answer.getIs_success())
            return answer;
        
        res = getCountNotificationLog(item);
        if(!res.getIs_success())
            return res;
        
        answer.setRes("[" + answer.getRes() + "," + res.getRes().substring(1, res.getRes().length() - 1) + "]");
        answer.setIs_success(res.getIs_success());
        
        return answer;
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
