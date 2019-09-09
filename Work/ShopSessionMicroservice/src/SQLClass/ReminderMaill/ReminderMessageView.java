package SQLClass.ReminderMaill;

import MyInterface.ModelInterface;
import java.sql.Timestamp;

public class ReminderMessageView implements ModelInterface {
    int id = -1; 
    String message_topic = "";
    Timestamp date_send; 
    
    public ReminderMessageView(int id, String message_topic, Timestamp date_send) {
        this.id = id;
        this.message_topic = message_topic;
        this.date_send = date_send;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), message_topic, date_send.toString()};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            "'" + message_topic.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            date_send.toString()
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "message_topic", "date_send"}; 
    }

    @Override
    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public String getMessage_topic() {
        return message_topic;
    }

    public void setMessage_topic(String message_topic) {
        this.message_topic = message_topic;
    }

    public Timestamp getDate_send() {
        return date_send;
    }

    public void setDate_send(Timestamp date_send) {
        this.date_send = date_send;
    }
}
