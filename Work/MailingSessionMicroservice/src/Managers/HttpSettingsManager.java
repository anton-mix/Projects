package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.GlobalVariables;
import SQLClass.HttpSettings.*;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Connection;
import java.sql.SQLException;

public class HttpSettingsManager implements ManagerInterface {

    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;
    String name_class = "HttpSettingsManager";
    
    public HttpSettingsManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.log_manager = log_manager;
        this.connection_with_db = connection_with_db;
        this.global_var = global_var;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(HttpSettings.class, new HttpSettingsConverter());
        builder.registerTypeAdapter(SocialIcons.class, new SocialIconsConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        builder.registerTypeAdapter(UserEmails.class, new UserEmailsConverter());
        gson = builder.create();
    }
    
    public ResData getHttpSettings(HttpSettings item) throws MyNotAuthorizedException {
        try {
            HttpSettings[] buff = gson.fromJson(Select("SELECT * FROM HttpSettings " + SqlFilter(item, item.getLimit(), item.getOffset()) + ";", connection_with_db, gson), HttpSettings[].class);
            for (int i = 0; i < buff.length;i++) {
                SocialIcons[] socialIcons = gson.fromJson(Select("SELECT * FROM SocialIcons WHERE http_settings_id = " + buff[i].getId() + " AND user_id = " + buff[i].getUser_id() + ";", connection_with_db, gson), SocialIcons[].class);
                for (int j = 0; j < socialIcons.length; j++) {
                    buff[i].getSocial_icons().add(socialIcons[j]);
                }
                UserEmails[] mails = gson.fromJson(Select("SELECT * FROM UserEmails WHERE http_settings_id = " + buff[i].getId() + " AND user_id = " + buff[i].getUser_id() + ";", connection_with_db, gson), UserEmails[].class);
                for (int j = 0; j < mails.length; j++) {
                    buff[i].getUser_emails().add(mails[j]);
                }
            }
            return new ResData(true, gson.toJson(buff, HttpSettings[].class));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getHttpSettings. Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getHttpSettings. Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getHttpSettings. Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getHttpSettings. Exception", name_class);
        }
    }
    
    public ResData AddHttpSettings(HttpSettings item) throws MyNotAuthorizedException {  
        try {
            ResData resReqwesr = getHttpSettings(new HttpSettings(item.getUser_id()));
            HttpSettings[] httpSettings = gson.fromJson(resReqwesr.getRes(), HttpSettings[].class);
            if(httpSettings.length >= 1) {
                return UpdateHttpSettings(item);
            } 
        }
        catch(Exception ex) {
            log_manager.CallError(ex.toString(), "AddHttpSettings. Error update defore add", name_class);
        }
        
        try {
            
            
            String res = Add(item, "HttpSettings", connection_with_db, gson, true);
            int http_settings_id = gson.fromJson(res, SQLResult.class).getId();
            
            for (int i = 0; i < item.getSocial_icons().size(); i++)  {
                item.getSocial_icons().get(i).setHttp_settings_id(http_settings_id);
                Add(item.getSocial_icons().get(i), "SocialIcons", connection_with_db, gson, true);
            }
            for (int i = 0; i < item.getUser_emails().size(); i++){
                item.getUser_emails().get(i).setHttp_settings_id(http_settings_id);
                Add(item.getUser_emails().get(i), "UserEmails", connection_with_db, gson, true);
            }
            return new ResData(true, res);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddHttpSettings. Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddHttpSettings. Exception in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddHttpSettings. Exception", name_class);
        }
    }
    
    public ResData AddUserEmails(UserEmails item) throws MyNotAuthorizedException {
        try {
            String res = Add(item, "UserEmails", connection_with_db, gson, true);
            return new ResData(true, res);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddHttpSettings. Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddHttpSettings. Exception in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddHttpSettings. Exception", name_class);
        }
    }
    
    public ResData DeleteHttpSettings(HttpSettings item) throws MyNotAuthorizedException {
        try {
            return new ResData(true,  Delete(item, "HttpSettings", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteHttpSettings. Sql add", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteHttpSettings. Sql add exception", name_class);
        }
    }
    
    public ResData UpdateHttpSettings(HttpSettings item) throws MyNotAuthorizedException {
        try {
            String res = ExecuteQuery("DELETE FROM SocialIcons where http_settings_id = " + item.getId() + ";", connection_with_db, gson);
            res = ExecuteQuery("DELETE FROM UserEmails where http_settings_id = " + item.getId() + ";", connection_with_db, gson);
            res = Update(item, "HttpSettings", connection_with_db, gson);
            for (int i = 0; i < item.getSocial_icons().size(); i++){
                item.getSocial_icons().get(i).setHttp_settings_id(item.getId());
                item.getSocial_icons().get(i).setUser_id(item.getUser_id());
                Add(item.getSocial_icons().get(i), "SocialIcons", connection_with_db, gson, true);
            }
            for (int i = 0; i < item.getUser_emails().size(); i++){
                item.getUser_emails().get(i).setHttp_settings_id(item.getId());
                item.getUser_emails().get(i).setUser_id(item.getUser_id());
                Add(item.getUser_emails().get(i), "UserEmails", connection_with_db, gson, true);
            }
            return new ResData(true, res);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateHttpSettings. Sql add", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateHttpSettings. Sql add exception", name_class);
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
