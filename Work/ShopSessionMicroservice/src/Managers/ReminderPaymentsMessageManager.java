package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.GlobalVariables;
import SQLClass.ReminderMaill.MessageImageLinks;
import SQLClass.ReminderMaill.MessageLinks;
import SQLClass.ReminderMaill.ReminderMessage;
import SQLClass.ReminderMaill.ReminderMessageConverter;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import SQLClass.User.UserSystemInfo;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Connection;
import java.sql.SQLException;

public class ReminderPaymentsMessageManager implements ManagerInterface {
    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;

    String name_class = "ReminderPaymentsMessageManager";

    public ReminderPaymentsMessageManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.global_var = global_var;
        this.connection_with_db = connection_with_db;
        this.log_manager = log_manager;
      
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(ReminderMessage.class, new ReminderMessageConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        gson = builder.create();      
    }
    
    public ResData getReminderPaymentsMessage(ReminderMessage item) throws MyNotAuthorizedException {
        try {
            ReminderMessage[] message = gson.fromJson(Select("SELECT * FROM ReminderPaymentsMailsTemplate " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson), ReminderMessage[].class);
            for(int i = 0; i < message.length; i++) {
                MessageLinks[] message_links = gson.fromJson(Select("SELECT link_address, link_name FROM ReminderPaymentsMailsLinks where user_id=" + message[i].getUser_id() + " AND message_id=" + message[i].getId(), connection_with_db, gson), MessageLinks[].class);
                MessageImageLinks[] image_links = gson.fromJson(Select("SELECT image_address, image_name FROM ReminderPaymentsMailsImageLinks where user_id=" + message[i].getUser_id() + " AND message_id=" + message[i].getId(), connection_with_db, gson), MessageImageLinks[].class);
                for(int j = 0; j < message_links.length; j++) 
                    message[i].getMessage_links().add(new MessageLinks(message_links[j].getLink_address(), message_links[j].getLink_name(), message_links[j].getLink_text_link()));
                for(int j = 0; j < image_links.length; j++)
                    message[i].getImage_links().add(new MessageImageLinks(image_links[j].getImage_address(), image_links[j].getImage_name(), message_links[j].getLink_text_link()));
            }
            return new ResData(true, gson.toJson(message));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getReminderPaymentsMessage(ReminderMessage). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getReminderPaymentsMessage(ReminderMessage). Exeption in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getReminderPaymentsMessage(ReminderMessage). Exception", name_class);
        }
    }
    
    public ResData AddReminderPaymentsMessage(ReminderMessage item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Add(item, "ReminderPaymentsMailsTemplate", connection_with_db, gson, true));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddReminderPaymentsMessage(ReminderMessage). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddReminderPaymentsMessage(ReminderMessage). Exeption in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddReminderPaymentsMessage(ReminderMessage). Exception", name_class);
        }
    }
    
    public ResData DeleteReminderPaymentsMessage(ReminderMessage item) throws MyNotAuthorizedException {
        ResData res;
        try {
            UserSystemInfo buff = gson.fromJson(getUserSystemInfo(new UserSystemInfo(item.getUser_id())).getRes(), UserSystemInfo[].class)[0];
            if(buff.getReminder_payment_id_1() != item.getReminder_id()&& buff.getReminder_payment_id_2() != item.getReminder_id() && buff.getReminder_payment_id_3() != item.getReminder_id()) {
                res = new ResData(true, Delete(item, "ReminderPaymentsMailsTemplate", connection_with_db, gson));
            }
            else {
                res = new ResData(true, gson.toJson(new SQLResult("Ok"), SQLResult.class));
            }
            return res;
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteReminderPaymentsMessage(ReminderMessage). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteReminderPaymentsMessage(ReminderMessage). Exception", name_class);
        }
    }
    
    public ResData UpdateReminderPaymentsMessage(ReminderMessage item) throws MyNotAuthorizedException{
        try {
            return new ResData(true, Update(item, "ReminderPaymentsMailsTemplate", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateReminderPaymentsMessage(ReminderMessage). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateReminderPaymentsMessage(ReminderMessage). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateReminderPaymentsMessage(ReminderMessage). Exception", name_class);
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
