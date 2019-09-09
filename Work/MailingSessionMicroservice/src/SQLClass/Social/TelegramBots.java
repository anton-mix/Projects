package SQLClass.Social;

import MyInterface.ModelInterface;

public class TelegramBots implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int type_id = -1;
    String bot_token = "";
    String bot_name = "";
    int enabled = 0;
    int banned = 0;
    
    public TelegramBots(int user_id, int type_id, String bot_token, String bot_name, int enabled, int banned) {
        this.user_id = user_id;
        this.type_id = type_id;
        this.bot_token = bot_token;
        this.bot_name = bot_name;
        this.enabled = enabled;
        this.banned = banned;
    }
    
    public TelegramBots(int id, int user_id, int type_id, String bot_token, String bot_name, int enabled, int banned) {
        this.id = id;
        this.user_id = user_id;
        this.type_id = type_id;
        this.bot_token = bot_token;
        this.bot_name = bot_name;
        this.enabled = enabled;
        this.banned = banned;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(type_id), bot_token, bot_name, String.valueOf(enabled), String.valueOf(banned)};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id), 
            String.valueOf(type_id), 
            "'" + bot_token.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + bot_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            String.valueOf(enabled),
            String.valueOf(banned)
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "type_id", "bot_token", "bot_name", "enabled", "banned"}; 
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

    public String getBot_token() {
        return bot_token;
    }

    public void setBot_token(String bot_token) {
        this.bot_token = bot_token;
    }

    public String getBot_name() {
        return bot_name;
    }

    public void setBot_name(String bot_name) {
        this.bot_name = bot_name;
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
}
