
package SQLClass.HttpSettings;

import MyInterface.ModelInterface;


public class SocialIcons implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int http_settings_id = -1;
    String icon_path = "";
    String link = "";
    
    public SocialIcons(int id, int user_id, int http_settings_id, String icon_path, String link){
       this.id = id;
       this.user_id = user_id;
       this.http_settings_id = http_settings_id;
       this.icon_path = icon_path;
       this.link = link;
    }
    
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(http_settings_id), icon_path, link};
    }

    @Override
    public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id),
            String.valueOf(http_settings_id), 
            "'" + icon_path.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + link.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
        };
    }

    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "http_settings_id", "icon_path", "link"}; 
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

    public int getHttp_settings_id() {
        return http_settings_id;
    }

    public void setHttp_settings_id(int http_settings_id) {
        this.http_settings_id = http_settings_id;
    }

    public String getLink() {
        return link;
    }

    public void setLink(String link) {
        this.link = link;
    }

    public String getIcon_path() {
        return icon_path;
    }

    public void setIcon_path(String icon_path) {
        this.icon_path = icon_path;
    } 
}
