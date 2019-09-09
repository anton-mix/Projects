package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.GlobalVariables;
import SQLClass.Product.ShortProduct;
import SQLClass.ReminderMaill.ReminderMessageView;
import SQLClass.ReminderMaill.ReminderMessageViewConverter;
import SQLClass.ReminderPayments.ReminderPaymentsModel;
import SQLClass.ReminderSurcharge.ReminderSurchargeModel;
import SQLClass.ReminderSurcharge.ReminderSurchargeModelConverter;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.User.UserSystemInfo;
import SQLClass.User.UserSystemInfoConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Connection;
import java.sql.SQLException;

public class ReminderSurchargeManager implements ManagerInterface {
    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;

    String name_class = "ReminderSurchargeManager";

    public ReminderSurchargeManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.global_var = global_var;
        this.connection_with_db = connection_with_db;
        this.log_manager = log_manager;
      
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(ReminderSurchargeModel.class, new ReminderSurchargeModelConverter());
        builder.registerTypeAdapter(ReminderMessageView.class, new ReminderMessageViewConverter());
        builder.registerTypeAdapter(UserSystemInfo.class, new UserSystemInfoConverter());
        gson = builder.create();      
    }
    
    public ResData getReminderSurcharge(ReminderSurchargeModel item) throws MyNotAuthorizedException {
        try {
            ReminderSurchargeModel[] items = gson.fromJson(Select("SELECT * FROM ReminderSurcharges " + SqlFilter(item, item.getLimit(), item.getOffset()) + ";", connection_with_db, gson), ReminderSurchargeModel[].class);
            for(int i = 0; i < items.length; i++) {
                ShortProduct[] products = gson.fromJson(Select("SELECT id, value, product_name FROM Products WHERE user_id = " + item.getUser_id() + " AND reminder_surcharge_id =  " + items[i].getId() + ";", connection_with_db, gson), ShortProduct[].class);
                for(int j = 0; j < products.length; j++) {
                    items[i].getProducts().add(products[j]);
                }
                
                ReminderMessageView[] messages = gson.fromJson(Select("SELECT id, message_topic, date_send FROM ReminderSurchargesMailsTemplate WHERE user_id = " + item.getUser_id() + " AND reminder_id =  " + items[i].getId() + ";", connection_with_db, gson), ReminderMessageView[].class);
                for(int j = 0; j < messages.length; j++) {
                    items[i].getMessages().add(messages[j]);
                }
                items[i].setCount(messages.length);
            }
            return new ResData(true, gson.toJson(items, ReminderSurchargeModel[].class));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getReminderSurcharge(ReminderSurchargeModel). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getReminderSurcharge(ReminderSurchargeModel). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getReminderSurcharge(ReminderSurchargeModel). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getReminderSurcharge(ReminderSurchargeModel). Exception", name_class);
        }
    }
    
    public ResData getCountReminderSurcharge(ReminderSurchargeModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT COUNT(*) as result FROM ReminderSurcharges " + SqlFilter(item, item.getLimit(), item.getOffset()) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getCountReminderSurcharge(ReminderPaymentsModel). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getCountReminderSurcharge(ReminderPaymentsModel). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getCountReminderSurcharge(ReminderPaymentsModel). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getCountReminderSurcharge(ReminderPaymentsModel). Exception", name_class);
        }
    }
    
    public ResData getSelectWithCountReminderPayments(ReminderSurchargeModel item){
        try{
            ResData res, answer;
            
            answer = getReminderSurcharge(item);
            if(!answer.getIs_success()) {
                return answer;
            }
            res = getCountReminderSurcharge(item);
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

    public ResData AddReminderSurcharge(ReminderSurchargeModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Add(item, "ReminderSurcharges", connection_with_db, gson, true));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddReminderSurcharge(ReminderSurchargeModel). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddReminderSurcharge(ReminderSurchargeModel). Exeption in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddReminderSurcharge(ReminderSurchargeModel). Exception", name_class);
        }
    }
    
    public ResData DeleteReminderSurcharge(ReminderSurchargeModel item) throws MyNotAuthorizedException {
        ResData res;
        try {
            UserSystemInfo buff = gson.fromJson(getUserSystemInfo(new UserSystemInfo(item.getUser_id())).getRes(), UserSystemInfo[].class)[0];
            if(buff.getReminder_surcharge_id_1()!= item.getId() && buff.getReminder_surcharge_id_2() != item.getId() && buff.getReminder_surcharge_id_3() != item.getId()) {
                res = new ResData(true, Delete(item, "ReminderSurcharges", connection_with_db, gson));
            }
            else {
                res =  new ResData(true, gson.toJson(new SQLResult("Ok"), SQLResult.class));
            }
            
            return res;
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteReminderSurcharge(ReminderSurchargeModel). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteReminderSurcharge(ReminderSurchargeModel). Exception", name_class);
        }
    }
    
    public ResData UpdateReminderSurcharge(ReminderSurchargeModel item) throws MyNotAuthorizedException{
        try {
            return new ResData(true, Update(item, "ReminderSurcharges", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateReminderSurcharge(ReminderSurchargeModel). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateReminderSurcharge(ReminderSurchargeModel). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateReminderSurcharge(ReminderSurchargeModel). Exception", name_class);
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
