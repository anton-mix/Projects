package SQLClass.Social;

import MyInterface.ModelInterface;

public class VkGroups implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int type_id = -1;
    String group_name = "";
    String group_token = "";
    int group_id = -1;
    int enabled = 0;
    int banned = 0;
   
    public VkGroups(int user_id, int type_id, String group_name, String group_token, int group_id, int enabled, int banned) {
        this.user_id = user_id;
        this.type_id = type_id;
        this.group_name = group_name;
        this.group_token = group_token;
        this.group_id = group_id;
        this.enabled = enabled;
        this.banned = banned;
    }
    
    public VkGroups(int id, int user_id, int type_id, String group_name, String group_token, int group_id, int enabled, int banned) {
        this.id = id;
        this.user_id = user_id;
        this.type_id = type_id;
        this.group_name = group_name;
        this.group_token = group_token;
        this.group_id = group_id;
        this.enabled = enabled;
        this.banned = banned;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(type_id), group_name, group_token, String.valueOf(group_id), String.valueOf(enabled), String.valueOf(banned)};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id), 
            String.valueOf(type_id), 
            "'" + group_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + group_token.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            String.valueOf(group_id),
            String.valueOf(enabled),
            String.valueOf(banned)
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "type_id", "group_name", "group_token", "group_id", "enabled", "banned"}; 
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

    public int getType_id() {
        return type_id;
    }

    public void setType_id(int type_id) {
        this.type_id = type_id;
    }

    public String getGroup_name() {
        return group_name;
    }

    public void setGroup_name(String group_name) {
        this.group_name = group_name;
    }

    public int getGroup_id() {
        return group_id;
    }

    public void setGroup_id(int group_id) {
        this.group_id = group_id;
    }

    public int getEnabled() {
        return enabled;
    }

    public void setEnabled(int enabled) {
        this.enabled = enabled;
    }

    public int getBanned() {
        return banned;
    }

    public void setBanned(int banned) {
        this.banned = banned;
    }

    public String getGroup_token() {
        return group_token;
    }

    public void setGroup_token(String group_token) {
        this.group_token = group_token;
    }
}
