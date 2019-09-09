
package SQLClass.Message;

import MyInterface.ModelInterface;
import java.sql.Timestamp;

public class ViewSingleMessage implements ModelInterface{
    int id = -1;
    int user_id = -1;
    String message_head = "";
    Timestamp date_send;
    String view_date = "";
    String view_time = "";
    String amount_receivers = "";
    String amount_opening = "";
    String amount_visiters = "";
    String sorted_by = "";
    
    public ViewSingleMessage(int id, int user_id, String message_head, Timestamp date_send, String view_date, String view_time,
        String amount_receivers, String amount_opening, String amount_visiters, String sorted_by){
        this.id = id;
        this.user_id = user_id;
        this.message_head = message_head;
        this.date_send = date_send;
        this.view_date = view_date;
        this.view_time = view_time;
        this.amount_receivers = amount_receivers;
        this.amount_opening = amount_opening;
        this.amount_visiters = amount_visiters;
        this.sorted_by = sorted_by;
    }
    
       
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), message_head, date_send.toString(),
            view_date, view_time};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id),
            String.valueOf(user_id), 
            "'" + message_head.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + date_send.toString().replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + view_date.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + date_send.toString().replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + view_date.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + view_time.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "message_topic", "date_send", "view_date", "view_time"}; 
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

    public String getMessage_head() {
        return message_head;
    }

    public void setMessage_head(String message_head) {
        this.message_head = message_head;
    }

    public Timestamp getDate_send() {
        return date_send;
    }

    public void setDate_send(Timestamp date_send) {
        this.date_send = date_send;
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

    public String getAmount_receivers() {
        return amount_receivers;
    }

    public void setAmount_receivers(String amount_receivers) {
        this.amount_receivers = amount_receivers;
    }

    public String getAmount_opening() {
        return amount_opening;
    }

    public void setAmount_opening(String amount_opening) {
        this.amount_opening = amount_opening;
    }

    public String getAmount_visiters() {
        return amount_visiters;
    }

    public void setAmount_visiters(String amount_visiters) {
        this.amount_visiters = amount_visiters;
    }

    public String getSorted_by() {
        return sorted_by;
    }

    public void setSorted_by(String sorted_by) {
        this.sorted_by = sorted_by;
    }
    
    
}
