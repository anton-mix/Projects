package QueueClass;

public class UserUpdatePassword {
    int user_id = -1;
    String old_password = "";
    String new_password = "";
    
    public UserUpdatePassword(int user_id, String old_password, String new_password) {
        this.user_id = user_id;
        this.old_password = old_password;
        this.new_password = new_password;
    }

    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public String getOld_password() {
        return old_password;
    }

    public void setOld_password(String old_password) {
        this.old_password = old_password;
    }

    public String getNew_password() {
        return new_password;
    }

    public void setNew_password(String new_password) {
        this.new_password = new_password;
    }
    
    
}
