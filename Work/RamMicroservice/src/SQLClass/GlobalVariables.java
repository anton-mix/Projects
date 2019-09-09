package SQLClass;

public class GlobalVariables {
    int port = -1;
    String queue_name = "", queue_name_log = "", host = "";
    public GlobalVariables(String host, String queue_name, String queue_name_log, int port) {
        this.host = host;
        this.queue_name = queue_name;
        this.queue_name_log = queue_name_log;
        this.port = port;
    }

    public String getQueue_name() {
        return queue_name;
    }

    public void setQueue_name(String queue_name) {
        this.queue_name = queue_name;
    }

    public String getQueue_name_log() {
        return queue_name_log;
    }

    public void setQueue_name_log(String queue_name_log) {
        this.queue_name_log = queue_name_log;
    }

    public String getHost() {
        return host;
    }

    public void setHost(String host) {
        this.host = host;
    }

    public int getPort() {
        return port;
    }

    public void setPort(int port) {
        this.port = port;
    }
}
