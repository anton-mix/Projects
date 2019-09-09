package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.GlobalVariables;
import SQLClass.Product.ShortProduct;
import SQLClass.ReminderMaill.ReminderMail;
import SQLClass.ReminderMaill.ReminderMailConverter;
import SQLClass.ReminderMaill.ReminderMessageView;
import SQLClass.ReminderMaill.ReminderMessageViewConverter;
import SQLClass.ReminderPayments.ReminderPaymentsModel;
import SQLClass.ReminderPayments.ReminderPaymentsModelConverter;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.User.UserSystemInfo;
import SQLClass.User.UserSystemInfoConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Connection;
import java.sql.SQLException;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.time.LocalDate;
import java.util.Date;

public class ReminderPaymentsManager implements ManagerInterface {
    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;

    String name_class = "ReminderPaymentsManager";

    public ReminderPaymentsManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.global_var = global_var;
        this.connection_with_db = connection_with_db;
        this.log_manager = log_manager;
      
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(ReminderPaymentsModel.class, new ReminderPaymentsModelConverter());
        builder.registerTypeAdapter(ReminderMessageView.class, new ReminderMessageViewConverter());
        builder.registerTypeAdapter(UserSystemInfo.class, new UserSystemInfoConverter());
        builder.registerTypeAdapter(ReminderMail.class, new ReminderMailConverter());
        gson = builder.create();      
    }
    
    public ResData getReminderPayments(ReminderPaymentsModel item) throws MyNotAuthorizedException {
        try {
            ReminderPaymentsModel[] items = gson.fromJson(Select("SELECT * FROM ReminderPayments " + SqlFilter(item, item.getLimit(), item.getOffset()) + ";", connection_with_db, gson), ReminderPaymentsModel[].class);
            for(int i = 0; i < items.length; i++) {
                ShortProduct[] products = gson.fromJson(Select("SELECT id, product_name FROM Products WHERE user_id = " + item.getUser_id() + " AND reminder_payment_id =  " + items[i].getId() + ";", connection_with_db, gson), ShortProduct[].class);
                for(int j = 0; j < products.length; j++) {
                    items[i].getProducts().add(products[j]);
                }
                
                ReminderMessageView[] messages = gson.fromJson(Select("SELECT id, message_topic, date_send FROM ReminderPaymentsMailsTemplate WHERE user_id = " + item.getUser_id() + " AND reminder_id =  " + items[i].getId() + ";", connection_with_db, gson), ReminderMessageView[].class);
                for(int j = 0; j < messages.length; j++) {
                    items[i].getMessages().add(messages[j]);
                }
                items[i].setCount(messages.length);
            }
            
            return new ResData(true, gson.toJson(items, ReminderPaymentsModel[].class));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getReminderPayments(ReminderPaymentsModel). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getReminderPayments(ReminderPaymentsModel). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getReminderPayments(ReminderPaymentsModel). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getReminderPayments(ReminderPaymentsModel). Exception", name_class);
        }
    }
    
    public ResData getCountReminderPayments(ReminderPaymentsModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT COUNT(*) as result FROM ReminderPayments " + SqlFilter(item, item.getLimit(), item.getOffset()) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getCountReminderPaymentsModel(ReminderPaymentsModel). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getCountReminderPaymentsModel(ReminderPaymentsModel). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getCountReminderPaymentsModel(ReminderPaymentsModel). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getCountReminderPaymentsModel(ReminderPaymentsModel). Exception", name_class);
        }
    }
    
    public ResData getSelectWithCountReminderPayments(ReminderPaymentsModel item){
        try{
            ResData res, answer;
            
            answer = getReminderPayments(item);
            if(!answer.getIs_success()) {
                return answer;
            }
            res = getCountReminderPayments(item);
            if(!res.getIs_success()) {
                return res;
            }
            
            answer.setRes("[" + answer.getRes() + "," + res.getRes().substring(1, res.getRes().length()-1)  + "]");
            answer.setIs_success(res.getIs_success());
            
            return answer;
        }
        catch(MyNotAuthorizedException ex){
            return log_manager.CallError(ex.toString(), "getSelectWithCountReminderPayments(ReminderPaymentsModel). MyNotAuthorizedException", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getSelectWithCountReminderPayments(ReminderPaymentsModel). Exeption", name_class);
        }
    }

    public ResData AddReminderPayments(ReminderPaymentsModel item) throws MyNotAuthorizedException {
        String result = "";
        int id = -1;
        try {
            result = Add(item, "ReminderPayments", connection_with_db, gson, true);
            id = gson.fromJson(result, SQLResult.class).getId();
            for (int i = 0; i < item.getProducts().size(); i++) {
                ExecuteQuery("UPDATE Products SET reminder_payment_id = " + id + " WHERE id = " + item.getProducts().get(i).getId() + ";", connection_with_db, gson);
            }
            
            return new ResData(true, Add(item, "ReminderPayments", connection_with_db, gson, true));
            
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddReminderPaymentsModel(ReminderPaymentsModel). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddReminderPaymentsModel(ReminderPaymentsModel). Exeption in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddReminderPaymentsModel(ReminderPaymentsModel). Exception", name_class);
        }
    }
    
    public ResData DeleteReminderPayments(ReminderPaymentsModel item) throws MyNotAuthorizedException {
        ResData res;
        try {
            UserSystemInfo buff = gson.fromJson(getUserSystemInfo(new UserSystemInfo(item.getUser_id())).getRes(), UserSystemInfo[].class)[0];
            if(buff.getReminder_payment_id_1() != item.getId() && buff.getReminder_payment_id_2() != item.getId() && buff.getReminder_payment_id_3() != item.getId()) {
                res = new ResData(true, Delete(item, "ReminderPayments", connection_with_db, gson));
            }
            else {
                res =  new ResData(true, gson.toJson(new SQLResult("Ok"), SQLResult.class));
            }
            return res;
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteReminderPayments(ReminderPaymentsModel). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteReminderPayments(ReminderPaymentsModel). Exception", name_class);
        }
    }
    
    public ResData UpdateReminderPayments(ReminderPaymentsModel item) throws MyNotAuthorizedException{
        ResData res = null;
        try {
            res = new ResData(true, Update(item, "ReminderPayments", connection_with_db, gson));
            for(int i =0; i < item.getProducts().size();i++) {
                ExecuteQuery("UPDATE Products SET reminder_payment_id = " + item.getId() + " WHERE id = " + item.getProducts().get(i).getId() + ";", connection_with_db, gson);
            }
            return res;
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateReminderPayments(ReminderPaymentsModel). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateReminderPayments(ReminderPaymentsModel). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateReminderPayments(ReminderPaymentsModel). Exception", name_class);
        }
    }
    
    public ResData addReminderPaymentsMail(ReminderMail item) {
        try {
            ReminderPaymentsModel buff = gson.fromJson(Select("SELECT * FROM ReminderPayments WHERE id = (select reminder_payment_id from Products where id = " + item.getId() + ");", connection_with_db, gson), ReminderPaymentsModel[].class)[0];
            
            for(int i = 0; i < buff.getMessages().size(); i++) {
                ReminderMail insert_item = new ReminderMail(item.getUser_id(), item.getContact_id(), buff.getMessages().get(i).getId(), item.getProduct_id(), Timestamp.valueOf("1970-01-02 0:0:0"));
                Add(insert_item, "ReminderPaymentsMails", connection_with_db, gson, true);
            }
            return new ResData(true, gson.toJson(new SQLResult("Ok"), SQLResult.class));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "addReminderPaymentsMail(ReminderMail). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "addReminderPaymentsMail(ReminderMail). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "addReminderPaymentsMail(ReminderMail). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "addReminderPaymentsMail(ReminderMail). Sql select Exception", name_class);
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
