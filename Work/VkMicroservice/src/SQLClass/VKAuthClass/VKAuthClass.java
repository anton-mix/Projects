package SQLClass.VKAuthClass;

public class VKAuthClass {
    int id = -1;
    int user_id = -1;
    String token = "";
    
    public VKAuthClass(int user_id, String token){
        this.id = id;
        this.token = token;
    }
    
    public VKAuthClass(int id, int user_id, String token){
        this.id = id;
        this.user_id = user_id;
        this.token = token;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public String getToken() {
        return token;
    }

    public void setToken(String token) {
        this.token = token;
    }

    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }
}
