package SQLClass;

public class GlobalVariables {
    int port = -1, db_port_mem = -1;
    String queue_name = "", queue_name_log = "", host = "";
    String db_ip_sql = "", db_name = "", db_user_name = "", db_password = "";
    String db_ip_mem = "", db_mem_user_name = "", db_mem_password = "";
    public GlobalVariables(String host, String queue_name, String queue_name_log, int port, String db_ip_sql, String db_name, 
            String db_user_name, String db_password, String db_ip_mem, int db_port_mem, String db_mem_user_name, String db_mem_password) {
        this.host = host;
        this.queue_name = queue_name;
        this.queue_name_log = queue_name_log;
        this.port = port;
        this.db_ip_sql = db_ip_sql;
        this.db_name = db_name;
        this.db_user_name = db_user_name;
        this.db_password = db_password;
        this.db_ip_mem = db_ip_mem;
        this.db_port_mem = db_port_mem;
        this.db_mem_user_name = db_mem_user_name;
        this.db_mem_password = db_mem_password;
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

    public String getDb_ip_sql() {
        return db_ip_sql;
    }

    public void setDb_ip_sql(String db_ip_sql) {
        this.db_ip_sql = db_ip_sql;
    }

    public String getDb_name() {
        return db_name;
    }

    public void setDb_name(String db_name) {
        this.db_name = db_name;
    }

    public String getDb_user_name() {
        return db_user_name;
    }

    public void setDb_user_name(String db_user_name) {
        this.db_user_name = db_user_name;
    }

    public String getDb_password() {
        return db_password;
    }

    public void setDb_password(String db_password) {
        this.db_password = db_password;
    }

    public int getPort() {
        return port;
    }

    public void setPort(int port) {
        this.port = port;
    }

    public int getDb_port_mem() {
        return db_port_mem;
    }

    public void setDb_port_mem(int db_port_mem) {
        this.db_port_mem = db_port_mem;
    }

    public String getDb_ip_mem() {
        return db_ip_mem;
    }

    public void setDb_ip_mem(String db_ip_mem) {
        this.db_ip_mem = db_ip_mem;
    }

    public String getDb_mem_user_name() {
        return db_mem_user_name;
    }

    public void setDb_mem_user_name(String db_mem_user_name) {
        this.db_mem_user_name = db_mem_user_name;
    }

    public String getDb_mem_password() {
        return db_mem_password;
    }

    public void setDb_mem_password(String db_mem_password) {
        this.db_mem_password = db_mem_password;
    }
    
}
