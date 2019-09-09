package SQLClass.User;

public class ShortUserEmail {
    int id = -1;
    String email = "";
    
    public ShortUserEmail() {
        id = -1;
        email = "";
    }
    
    public ShortUserEmail(int id, String email) {
        this.id = id;
        this.email = email;
    }

    public String getEmail() {
        return email;
    }

    public void setEmail(String email) {
        this.email = email;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }
}
