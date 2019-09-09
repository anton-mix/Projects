package QueueClass;

public class ResQueueItem {
    int id = -1;
    String session_key = "";
    
    public ResQueueItem(int id,String session_key) {
        this.id = id;
        this.session_key = session_key;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public String getSession_key() {
        return session_key;
    }

    public void setSession_key(String session_key) {
        this.session_key = session_key;
    }
    
    
}
