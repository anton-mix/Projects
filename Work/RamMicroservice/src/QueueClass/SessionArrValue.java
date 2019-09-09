package QueueClass;

import java.sql.Timestamp;

public class SessionArrValue {
    int id = -1;
    Timestamp date;
    
    public SessionArrValue(int id, Timestamp date) {
        this.id = id;
        this.date = date;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public Timestamp getDate() {
        return date;
    }

    public void setDate(Timestamp date) {
        this.date = date;
    }
}
