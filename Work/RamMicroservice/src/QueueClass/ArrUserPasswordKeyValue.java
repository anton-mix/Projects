package QueueClass;

import java.sql.Timestamp;

public class ArrUserPasswordKeyValue {
    Timestamp date;
    String login = "";
    
    public ArrUserPasswordKeyValue(Timestamp date, String login) {
        this.date = date;
        this.login = login;
    }

    public Timestamp getDate() {
        return date;
    }

    public void setDate(Timestamp date) {
        this.date = date;
    }

    public String getLogin() {
        return login;
    }

    public void setLogin(String login) {
        this.login = login;
    }
     
}
