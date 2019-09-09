package QueueClass;

import java.sql.Timestamp;

public class UserArrValue {
   Timestamp date;
   UserAuthModel user;
   
   public UserArrValue(Timestamp date, UserAuthModel user) {
       this.date = date;
       this.user = user;
   }

    public Timestamp getDate() {
        return date;
    }

    public void setDate(Timestamp date) {
        this.date = date;
    }

    public UserAuthModel getUser() {
        return user;
    }

    public void setUser(UserAuthModel user) {
        this.user = user;
    }
}
