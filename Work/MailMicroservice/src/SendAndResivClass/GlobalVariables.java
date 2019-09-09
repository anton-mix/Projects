package SendAndResivClass;

public class GlobalVariables {
    String queue_name = "", queue_name_log = "", host = "";
    String username = "", password = "", hostMail = "";
    int port = -1;
    public GlobalVariables(String host, String queue_name, String queue_name_log, int port, String username, String password, String hostMail) {
        this.host = host;
        this.queue_name = queue_name;
        this.queue_name_log = queue_name_log;
        this.port = port;
        this.username = username;
        this.password = password;
        this.hostMail = hostMail;
    }

    public String getQueueName() {
        return queue_name;
    }

    public void setQueueName(String queue_name) {
        this.queue_name = queue_name;
    }

    public String getQueueNameLog() {
        return queue_name_log;
    }

    public void setQueueNameLog(String queue_name_log) {
        this.queue_name_log = queue_name_log;
    }

    public String getHost() {
        return host;
    }

    public void setHost(String host) {
        this.host = host;
    }

    public String getUsername() {
        return username;
    }

    public void setUsername(String username){
        this.username = username;
    }

    public String getPassword() {
        return password;
    }

    public void setPassword(String password) {
        this.password = password;
    }

    public String getHostMail() {
        return hostMail;
    }

    public void setHostMail(String hostMail) {
        this.hostMail = hostMail;
    }

    public int getPort() {
        return port;
    }

    public void setPort(int port) {
        this.port = port;
    }
}
