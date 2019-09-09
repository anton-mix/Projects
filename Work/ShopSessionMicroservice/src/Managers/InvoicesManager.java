package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.Invoice.FilterInvoice;
import SQLClass.Invoice.InvoiceModel;
import SQLClass.GlobalVariables;
import SQLClass.Invoice.InvoiceModelConverter;
import SQLClass.Invoice.ViewInvoice;
import SQLClass.Invoice.ViewInvoiceConverter;
import SQLClass.Invoice.InvoiceToProductsItem;
import SQLClass.Product.ShortProduct;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import java.sql.Connection;
import java.sql.SQLException;

public class InvoicesManager implements ManagerInterface {
    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;
    
    String name_class = "InvoicesManager";
    
    public InvoicesManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.log_manager = log_manager;
        this.connection_with_db = connection_with_db;
        this.global_var = global_var;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(InvoiceModel.class, new InvoiceModelConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        builder.registerTypeAdapter(ViewInvoice.class, new ViewInvoiceConverter());
        gson = builder.create();
    }

    //Select
    
    public ResData getInvoices() {
        try {
            return new ResData(true, Select("SELECT * FROM Invoices;", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getInvoices(). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getInvoices(). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getInvoices(). Exception", name_class);
        }
    }

    public ResData getInvoice(ViewInvoice item) throws MyNotAuthorizedException {
        try {
            InvoiceModel[] invoice = gson.fromJson(Select("SELECT * FROM Invoices " + SqlFilter(item,0,0) + ";", connection_with_db, gson), InvoiceModel[].class);
            for(int i = 0; i < invoice.length; i++) {
                ShortProduct[] product = gson.fromJson(Select("SELECT c.id, c.value, c.product_name FROM Invoices a, Invoices_to_Products b, Products c WHERE b.user_id = " + invoice[i].getUser_id() + " AND a.id = " + invoice[i].getId() + " AND a.id = b.invoice_id AND b.product_id = c.id;", connection_with_db, gson), ShortProduct[].class);
                for(int j = 0; j < product.length; j++) {
                    invoice[i].getProducts().add(product[j]);
                }
            }
            return new ResData(true, gson.toJson(invoice));
            
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getInvoice(ViewInvoice). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getInvoice(ViewInvoice). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getInvoice(ViewInvoice). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getInvoices(ViewInvoice). Exceptions", name_class);
        }
    }
    
    public ResData getViewInvoices (FilterInvoice item) throws MyNotAuthorizedException {
        try {
            ViewInvoice[] invoice = gson.fromJson(Select("SELECT * FROM Invoices " + SqlFilter(item, item.getLimit(), item.getOffset()) + ";", connection_with_db, gson), ViewInvoice[].class);
            for(int i = 0; i < invoice.length; i++) {
                ShortProduct[] product = gson.fromJson(Select("SELECT c.id, c.value,c.product_name FROM Invoices a, Invoices_to_Products b, Products c WHERE b.user_id = " + invoice[i].getUser_id() + " AND a.id = " + invoice[i].getId() + " AND a.id = b.invoice_id AND b.product_id = c.id;", connection_with_db, gson), ShortProduct[].class);
                for(int j = 0; j < product.length; j++) {
                    invoice[i].getProducts().add(product[j]);
                }
            }
            return new ResData(true, gson.toJson(invoice));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getViewInvoices(FilterInvoice). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getViewInvoices(FilterInvoice). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getViewInvoices(FilterInvoice). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getViewInvoices(FilterInvoice). Exceptions", name_class);
        } 
    }

    public ResData getCountInvoices(FilterInvoice item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT COUNT(*) as result FROM Invoices " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getCountInvoices(FilterInvoice). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getCountInvoices(FilterInvoice). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getCountInvoices(FilterInvoice). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getCountInvoices(FilterInvoice). Exceptions", name_class);
        } 
    }

    //Add
    
    public ResData AddInvoice(InvoiceModel item) throws MyNotAuthorizedException {
        SQLResult res;
        String result;
        try {
            result = Add(item, "Invoices", connection_with_db, gson, true);
            res = gson.fromJson(result, SQLResult.class);
            for(int i = 0; i < item.getProducts().size(); i++) {
                InvoiceToProductsItem pitem = new InvoiceToProductsItem(res.getUser_id(),res.getId(),item.getProducts().get(i).getId());
                Add(pitem, "Invoices_to_Products", connection_with_db, gson, true);
            }
            return new ResData(true, result);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddInvoice(InvoiceModel). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddInvoice(InvoiceModel). Exeption in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddInvoice(InvoiceModel). Exceptions", name_class);
        }
    }
    
    //Delete
    
    public ResData DeleteInvoice(ViewInvoice item) throws MyNotAuthorizedException {
        try {
            ExecuteQuery("DELETE FROM InvoiceEvents WHERE user_id = " + item.getUser_id() + " AND invoice_id = " + item.getId() + ";", connection_with_db, gson);
            return new ResData(true, Delete(item, "Invoices", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteInvoice(ViewInvoice). Sql delete", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "DeleteInvoice(ViewInvoice). Sql all delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteInvoice(ViewInvoice). Exceptions", name_class);
        }
    }
    
    //Update
    
    public ResData UpdateInvoice(InvoiceModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Update(item, "Invoices", connection_with_db, gson));
  
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateInvoice(InvoiceModel). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateInvoice(InvoiceModel). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateInvoice(InvoiceModel). Exceptions", name_class);
        }
    }

    public ResData Invoice_select_with_count(FilterInvoice item){
        try {
            ResData res, answer;
            
            answer = getViewInvoices(item);
            if(!answer.getIs_success())
                return answer;
            
            res = getCountInvoices(item);
            if(!res.getIs_success())
                return res;
            
            answer.setRes("[" + answer.getRes() + "," + res.getRes().substring(1, res.getRes().length() - 1)  + "]");
            answer.setIs_success(res.getIs_success());
            
            return answer;
        }
        catch(MyNotAuthorizedException ex)  {
            return log_manager.CallError(ex.toString(), "Invoice_select_with_count(FilterInvoice). MyNotAuthorizedException", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "Invoice_select_with_count(FilterInvoice). Exceptions", name_class);
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
