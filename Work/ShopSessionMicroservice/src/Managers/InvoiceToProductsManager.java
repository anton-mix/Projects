package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.GlobalVariables;
import SQLClass.Invoice.InvoiceToProductsItem;
import SQLClass.Invoice.InvoiceToProductsItemConverter;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Connection;
import java.sql.SQLException;

public class InvoiceToProductsManager implements ManagerInterface {
    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;

    String name_class = "InvoiceEventManager";

    public InvoiceToProductsManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.global_var = global_var;
        this.connection_with_db = connection_with_db;
        this.log_manager = log_manager;
      
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(InvoiceToProductsItem.class, new InvoiceToProductsItemConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        gson = builder.create();      
    }
    
    public ResData getInvoiceToProductsItems(InvoiceToProductsItem item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT * FROM Invoices_to_Products " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getInvoiceToProductsItems(InvoiceToProductsItem). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getInvoiceToProductsItems(InvoiceToProductsItem). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getInvoiceToProductsItems(InvoiceToProductsItem). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getInvoiceToProductsItems(InvoiceToProductsItem). Exception", name_class);
        }
    }
    
    public ResData getCountInvoiceToProductsItems(InvoiceToProductsItem item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT COUNT(*) as result FROM Invoices_to_Products " + SqlFilter(item,0,0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getCountInvoiceToProductsItems(InvoiceToProductsItem). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getCountInvoiceToProductsItems(InvoiceToProductsItem). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getCountInvoiceToProductsItems(InvoiceToProductsItem). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getCountInvoiceToProductsItems(InvoiceToProductsItem). Exception", name_class);
        }
    }

    public ResData AddInvoiceToProductsItem(InvoiceToProductsItem item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Add(item, "Invoices_to_Products", connection_with_db, gson, true));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddInvoiceToProductsItem(InvoiceToProductsItem). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddInvoiceToProductsItem(InvoiceToProductsItem). Exeption in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddInvoiceToProductsItem(InvoiceToProductsItem). Exception", name_class);
        }
    }
    
    public ResData DeleteInvoiceToProductsItem(InvoiceToProductsItem item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Delete(item, "Invoices_to_Products", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteInvoiceToProductsItem(InvoiceToProductsItem). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddInvoiceToProductsItem(InvoiceToProductsItem). Exception", name_class);
        }
    }
    
    public ResData UpdateInvoiceToProductsItem(InvoiceToProductsItem item) throws MyNotAuthorizedException{
        try {
            return new ResData(true, Update(item, "Invoices_to_Products", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateInvoiceToProductsItem(InvoiceToProductsItem). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateInvoiceToProductsItem(InvoiceToProductsItem). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateInvoiceToProductsItem(InvoiceToProductsItem). Exception", name_class);
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
