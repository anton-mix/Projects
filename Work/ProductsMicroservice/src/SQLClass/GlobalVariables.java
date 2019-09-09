package SQLClass;

import java.util.ArrayList;
import java.util.List;

public class GlobalVariables {
    String db_ip_sql = "", db_name = "", db_user_name = "", db_password = "";
    String queue_name = "", queue_name_log = "", host = "";
    public GlobalVariables(String host, String queue_name, String queue_name_log, String db_ip_sql, String db_name, String db_user_name, String db_password) {
        this.host = host;
        this.queue_name = queue_name;
        this.queue_name_log = queue_name_log;
        this.db_ip_sql = db_ip_sql;
        this.db_name = db_name;
        this.db_user_name = db_user_name;
        this.db_password = db_password;
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
    
     public String getDbIpSql() {
        return db_ip_sql;
    }

    public void setDbIpSql(String db_ip_sql) {
        this.db_ip_sql = db_ip_sql;
    }

    public String getDbName() {
        return db_name;
    }

    public void setDbName(String db_name) {
        this.db_name = db_name;
    }

    public String getDbUserName() {
        return db_user_name;
    }

    public void setDbUserName(String db_user_name) {
        this.db_user_name = db_user_name;
    }

    public String getDbPassword() {
        return db_password;
    }

    public void setDbPassword(String db_password) {
        this.db_password = db_password;
    }
    
    public static List<String[]> Parse_json(String json) {
        String[] buff1, buff2;
        List<String[]> Answer = new ArrayList<String[]>();
        json = json.substring(2, json.length()-1);
        
        buff1 = json.split(",");

        for(int i = 0; i < buff1.length; i ++) {
            
            if("ul".equals(buff1[i])) {
                continue;
            }
            
            buff2 = buff1[i].split(":", 2);
            buff2[0] = buff2[0].substring(1,buff2[0].length()-1);
            
            if(buff2[1].indexOf("\"") != -1 && buff2[1].indexOf("\"") == 0) {
                buff2[1] = buff2[1].substring(1,buff2[1].length()-1);
            }

            Answer.add(new String[] {buff2[0],buff2[1]}); 
        }  
        return Answer;
    }
}
