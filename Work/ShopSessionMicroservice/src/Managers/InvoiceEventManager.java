package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.GlobalVariables;
import SQLClass.Invoice.InvoiceEvent;
import SQLClass.Invoice.InvoiceEventConverter;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Connection;
import java.sql.SQLException;

public class InvoiceEventManager implements ManagerInterface {
    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;
    
    String name_class = "InvoiceEventManager";
    
    public InvoiceEventManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.global_var = global_var;
        this.connection_with_db = connection_with_db;
        this.log_manager = log_manager;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(InvoiceEvent.class, new InvoiceEventConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        gson = builder.create();
    }
    
    public ResData getInvoiceEvents (InvoiceEvent item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT * FROM InvoiceEvents " + SqlFilter(item,0,0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getEventInvoice(InvoiceEvent). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getEventInvoice(InvoiceEvent). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getEventInvoice(InvoiceEvent). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getEventInvoice(InvoiceEvent). Exception", name_class);
        }
    }

    public ResData getCountInvoiceEvents(InvoiceEvent item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT COUNT(*) as result FROM InvoiceEvents " + SqlFilter(item,0,0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getCountInvoiceEvents(InvoiceEvent). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getCountInvoiceEvents(InvoiceEvent). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getCountInvoiceEvents(InvoiceEvent). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getCountInvoiceEvents(InvoiceEvent). Exception", name_class);
        }
    }
    

    public ResData AddInvoiceEvents(InvoiceEvent item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Add(item, "InvoiceEvents", connection_with_db, gson, true));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddInvoiceEvents(InvoiceEvent). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddInvoiceEvents(InvoiceEvent). Exeption in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddInvoiceEvents(InvoiceEvent). Exception", name_class);
        }
    }
    
    public ResData DeleteInvoiceEvents(InvoiceEvent item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Delete(item, "InvoiceEvents", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteInvoiceEvents(InvoiceEvent). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteInvoiceEvents(InvoiceEvent). Exception", name_class);
        }
    }
    
    public ResData UpdateInvoiceEvents(InvoiceEvent item) throws MyNotAuthorizedException{
        try {
            return new ResData(true, Update(item, "InvoiceEvents", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateInvoiceEvents(InvoiceEvent). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateInvoiceEvents(InvoiceEvent). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateInvoiceEvents(InvoiceEvent). Exception", name_class);
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
