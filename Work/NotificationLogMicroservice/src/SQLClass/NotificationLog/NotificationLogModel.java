package SQLClass.NotificationLog;

import MyInterface.ModelInterface;
import java.sql.Timestamp;

public class NotificationLogModel implements ModelInterface {
    int id = -1;
    int user_id = -1;
    Timestamp date;
    String content = "";
    
    public NotificationLogModel(int id, int user_id, Timestamp date, String content) {
        this.id = id;
        this.user_id = user_id;
        this.date = date;
        this.content = content;
    }
    
    public NotificationLogModel(int user_id, Timestamp date, String content) {
        this.user_id = user_id;
        this.date = date;
        this.content = content;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), date.toString(), content};
    }
    
    @Override
    public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id),
            String.valueOf(user_id), 
            "'" + date.toString() + "'",
            "'" + content.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
            };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "date", "content"}; 
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

    public Timestamp getDate() {
        return date;
    }

    public void setDate(Timestamp date) {
        this.date = date;
    }

    public String getContent() {
        return content;
    }

    public void setContent(String content) {
        this.content = content;
    }
    
    
}
