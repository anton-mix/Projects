package RestObject;
 
public class GlobalVariables {
    String host = "", queue_name = "";
    int port = -1;
    public GlobalVariables(String host, String queue_name, int port) {
        this.host = host;
        this.port = port;
        this.queue_name = queue_name;
    }

    public String getHost() {
        return host;
    }

    public void setHost(String host) {
        this.host = host;
    }

    public String getQueue_name() {
        return queue_name;
    }

    public void setQueue_name(String queue_name) {
        this.queue_name = queue_name;
    }

    public int getPort() {
        return port;
    }

    public void setPort(int port) {
        this.port = port;
    }
}
