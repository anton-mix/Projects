
package SQLClass.Message;

import MyInterface.ModelInterface;
import java.sql.Timestamp;


public class FilterMessage implements ModelInterface{
    int user_id = -1;
    int message_num = -1;
    Timestamp date_send_begin;
    Timestamp date_send_end; 
    String message_topic = "";
    int limit = 0;
    int offset = 0; 
    
    public FilterMessage(int user_id, int message_num, Timestamp date_send_begin, Timestamp date_send_end,
            String message_topic, int limit, int offset) {
        this.user_id = user_id;
        this.message_num = message_num;
        this.date_send_begin = date_send_begin;
        this.date_send_end = date_send_end;
        this.message_topic = message_topic;
        this.limit = limit;
        this.offset = offset;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(user_id), message_topic, String.valueOf(message_num), date_send_begin.toString(),
            date_send_end.toString(), String.valueOf(limit), String.valueOf(offset)};
    }
    
    @Override
    public String[] getSqlArray() {
        return new String[] {
            String.valueOf(user_id), 
            "'" + message_topic.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"user_id", "message_topic" ,"message_num", "date_send_begin", "date_send_end", "limit", "offset"}; 
    }

    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public int getMessage_num() {
        return message_num;
    }

    public void setMessage_num(int message_num) {
        this.message_num = message_num;
    }

    public Timestamp getDate_send_begin() {
        return date_send_begin;
    }

    public void setDate_send_begin(Timestamp date_send_begin) {
        this.date_send_begin = date_send_begin;
    }

    public Timestamp getDate_send_end() {
        return date_send_end;
    }

    public void setDate_send_end(Timestamp date_send_end) {
        this.date_send_end = date_send_end;
    }

    public String getMessage_topic() {
        return message_topic;
    }

    public void setMessage_topic(String message_topic) {
        this.message_topic = message_topic;
    }

    @Override
    public int getLimit() {
        return limit;
    }

    public void setLimit(int limit) {
        this.limit = limit;
    }

    @Override
    public int getOffset() {
        return offset;
    }

    public void setOffset(int offset) {
        this.offset = offset;
    }
    
    
}
