package SQLClass;

public class SQLResult {
    String result = "";
    int id = -1;
    int user_id = -1;
    
    public SQLResult(String result, int id, int user_id) {
        this.result = result;
        this.id = id;
        this.user_id = user_id;
    }
    
    public SQLResult(String result) {
        this.result = result;
        this.id = -1;
        this.user_id = -1;
    }

    public String getResult() {
        return result;
    }

    public void setResult(String result) {
        this.result = result;
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
    
    
}
