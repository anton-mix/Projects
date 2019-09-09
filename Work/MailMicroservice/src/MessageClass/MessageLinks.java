
package MessageClass;

import MyInterface.ModelInterface;

public class MessageLinks implements ModelInterface{
    int id = -1;
    int user_id = -1;
    int message_id = -1;
    String link_address = "";
    String link_name = "";
    
    public MessageLinks(int id, int user_id, int message_id, String link_address, String link_name) {
        this.id = id;
        this.user_id = user_id;
        this.message_id = message_id;
        this.link_address = link_address;
        this.link_name = link_name;
    }
    
    public MessageLinks(int user_id, int message_id, String link_address, String link_name) {
        this.user_id = user_id;
        this.message_id = message_id;
        this.link_address = link_address;
        this.link_name = link_name;
    }
    
    public MessageLinks(String link_address, String link_name) {
        this.link_address = link_address;
        this.link_name = link_name;
    }
        
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(message_id), 
            link_address, link_name};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id),
            String.valueOf(message_id), 
            "'" + link_address.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + link_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "message_id", "link_address", "link_name"}; 
    }

    @Override
    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public int getMessage_id() {
        return message_id;
    }

    public void setMessage_id(int message_id) {
        this.message_id = message_id;
    }

    public String getLink_address() {
        return link_address;
    }

    public void setLink_address(String link_address) {
        this.link_address = link_address;
    }

    public String getLink_name() {
        return link_name;
    }

    public void setLink_name(String link_name) {
        this.link_name = link_name;
    }
    
    
    
}
