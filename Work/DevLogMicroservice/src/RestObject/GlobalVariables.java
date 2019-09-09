package RestObject;
 
public class GlobalVariables {
    String host = "", queue_name = "", Queue_name_log = "";
    int port = -1;
    public GlobalVariables(String host, String queue_name, String Queue_name_log, int port) {
        this.host = host;
        this.port = port;
        this.queue_name = queue_name;
        this.Queue_name_log = Queue_name_log;
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

    public String getQueue_name_log() {
        return Queue_name_log;
    }

    public void setQueue_name_log(String Queue_name_log) {
        this.Queue_name_log = Queue_name_log;
    }

    public int getPort() {
        return port;
    }

    public void setPort(int port) {
        this.port = port;
    }
}
