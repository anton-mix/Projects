package RestObject;

import java.sql.Timestamp;

public class DevLog {
    int id = -1;
    int user_id = -1;
    Timestamp date;
    String theme = "";
    String content = "";
    
    public DevLog(int id, int user_id, Timestamp date, String theme, String content) {
        this.id = id;
        this.user_id = user_id;
        this.date = date;
        this.theme = theme;
        this.content = content;
    }
    
    public DevLog(int user_id, Timestamp date, String theme, String content) {
        this.user_id = user_id;
        this.date = date;
        this.theme = theme;
        this.content = content;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

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

    public String getTheme() {
        return theme;
    }

    public void setTheme(String theme) {
        this.theme = theme;
    }

    public String getContent() {
        return content;
    }

    public void setContent(String content) {
        this.content = content;
    }
    
    @Override
    public String toString() {
        return "user_id: " + user_id + " date: " + date.toString() + " theme: " + theme + " content: " + content;
    }
}
