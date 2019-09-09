package SQLClass.Invoice;

import MyInterface.ModelInterface;
import java.sql.Timestamp;

public class InvoiceEvent implements ModelInterface {
    int id = -1; 
    int user_id = -1; 
    int invoice_id = -1; 
    int type = -1; 
    String is_system = "";
    String title = ""; 
    Timestamp date;
    String owner = ""; 
    String is_call_back = "";
    String call_back = ""; 
    String comment = "";
    String view_date = "";
    String view_time = "";
    
    public InvoiceEvent(int id, int user_id, int invoice_id, int type, String is_system, String title, Timestamp date, String owner, String is_call_back, String call_back, String comment, String view_date, String view_time){
        this.id = id; 
        this.user_id = user_id; 
        this.invoice_id = invoice_id; 
        this.type = type; 
        this.is_system = is_system;
        this.title = title; 
        this.date = date; 
        this.owner = owner;
        this.is_call_back = is_call_back;
        this.call_back = call_back; 
        this.comment = comment;
        this.view_date = view_date;
        this.view_time = view_time;
    }
    
    public InvoiceEvent(int user_id, int invoice_id, int type, String is_system, String title, Timestamp date, String owner, String is_call_back, String call_back, String comment) {
        
        this.user_id = user_id; 
        this.invoice_id = invoice_id; 
        this.type = type; 
        this.is_system = is_system;
        this.title = title; 
        this.date = date; 
        this.owner = owner;
        this.is_call_back = is_call_back;
        this.call_back = call_back; 
        this.comment = comment;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(invoice_id), String.valueOf(type), is_system, title, date.toString(), 
            owner, is_call_back, call_back, comment, view_date, view_time};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id), 
            String.valueOf(invoice_id), 
            String.valueOf(type), 
            "'" + is_system.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + title.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + date.toString() + "'", 
            "'" + owner.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + is_call_back.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + call_back.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + comment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
        }; 
     }
     
    @Override
    public String[] getKeyArray() {
        return new String[] { "id", "user_id", "invoice_id", "type", "is_system","title", "date", "owner", "is_call_back", "call_back", "comment",
            "view_date", "view_time" }; 
    }

    public String getIs_system() {
        return is_system;
    }

    public void setIs_system(String is_system) {
        this.is_system = is_system;
    }
    
    @Override
    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public int getInvoice_id() {
        return invoice_id;
    }

    public void setInvoice_id(int invoice_id) {
        this.invoice_id = invoice_id;
    }

    public int getType() {
        return type;
    }

    public void setType(int type) {
        this.type = type;
    }

    public String getTitle() {
        return title;
    }

    public void setTitle(String title) {
        this.title = title;
    }

    public Timestamp getDate() {
        return date;
    }

    public void setDate(Timestamp date) {
        this.date = date;
    }

    public String getOwner() {
        return owner;
    }

    public void setOwner(String owner) {
        this.owner = owner;
    }

    public String getIs_call_back() {
        return is_call_back;
    }

    public void setIs_call_back(String is_call_back) {
        this.is_call_back = is_call_back;
    }
    
    public String getCall_back() {
        return call_back;
    }

    public void setCall_back(String call_back) {
        this.call_back = call_back;
    }

    public String getComment() {
        return comment;
    }

    public void setComment(String comment) {
        this.comment = comment;
    }

    public String getView_date() {
        return view_date;
    }

    public void setView_date(String view_date) {
        this.view_date = view_date;
    }

    public String getView_time() {
        return view_time;
    }

    public void setView_time(String view_time) {
        this.view_time = view_time;
    }
    
}
