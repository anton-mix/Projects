package SQLClass;

public class GlobalVariables {
    String dbIpSql, dbName, dbUserName, dbPassword;
    String queueNameLog, host;
    int port;
    String windows_base_path_image = "g://desctop//Project 3//Project 35//Files//Images//";
    String unix_base_path_image = "//home//server//Files//Images//";
    public GlobalVariables() {
        host = "141.8.198.38";
        queueNameLog = "log_queue";
        port = 5672;
        dbIpSql = "141.8.198.38";
        dbName = "UTP1";
        dbUserName = "utp_1_db";
        dbPassword = "db23db23";
    }
    
     public String getDbIpSql() {
        return dbIpSql;
     }

    public String getQueueNameLog() {
        return queueNameLog;
    }

    public void setQueueNameLog(String queueNameLog) {
        this.queueNameLog = queueNameLog;
    }

    public String getHost() {
        return host;
    }

    public void setHost(String host) {
        this.host = host;
    }

    public void setDbIpSql(String dbIpSql) {
        this.dbIpSql = dbIpSql;
    }

    public String getDbName() {
        return dbName;
    }

    public void setDbName(String dbName) {
        this.dbName = dbName;
    }

    public String getDbUserName() {
        return dbUserName;
    }

    public void setDbUserName(String dbUserName) {
        this.dbUserName = dbUserName;
    }

    public String getDbPassword() {
        return dbPassword;
    }

    public void setDbPassword(String dbPassword) {
        this.dbPassword = dbPassword;
    }

    public int getPort() {
        return port;
    }

    public void setPort(int port) {
        this.port = port;
    }

    public String getWindows_base_path_image() {
        return windows_base_path_image;
    }

    public void setWindows_base_path_image(String windows_base_path_image) {
        this.windows_base_path_image = windows_base_path_image;
    }

    public String getUnix_base_path_image() {
        return unix_base_path_image;
    }

    public void setUnix_base_path_image(String unix_base_path_image) {
        this.unix_base_path_image = unix_base_path_image;
    }

}
