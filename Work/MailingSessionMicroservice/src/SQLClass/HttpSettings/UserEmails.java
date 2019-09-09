
package SQLClass.HttpSettings;

import MyInterface.ModelInterface;


public class UserEmails implements ModelInterface {
    int id = -1; 
    int user_id = -1; 
    int http_settings_id = -1;
    String sender_name = "";
    String sender_email = "";
    String mailing_name = "";
    String main_url = "";
    
    public UserEmails(){
        this.id = -1; 
        this.user_id = -1; 
        this.http_settings_id = -1;
        this.sender_name = "";
        this.sender_email = "";
        this.mailing_name = "";
        this.main_url = "";
    }
    
    public UserEmails(int id, int user_id, int http_settings_id, String sender_name, String sender_email,
            String mailing_name, String main_url){
       this.id = id;
       this.user_id = user_id;
       this.http_settings_id = http_settings_id;
       this.sender_name = sender_name;
       this.sender_email = sender_email;
       this.mailing_name = mailing_name;
       this.main_url = main_url;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(http_settings_id), sender_name, 
            sender_email, mailing_name, main_url};
    }

    @Override
    public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id),
            String.valueOf(http_settings_id), 
            "'" + sender_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + sender_email.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + mailing_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + main_url.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
        };
    }

    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "http_settings_id", "sender_name", "sender_email", "mailing_name", "main_url"}; 
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

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

    public String getSender_name() {
        return sender_name;
    }

    public void setSender_name(String sender_name) {
        this.sender_name = sender_name;
    }

    public String getSender_email() {
        return sender_email;
    }

    public void setSender_email(String sender_email) {
        this.sender_email = sender_email;
    }

    public String getMailing_name() {
        return mailing_name;
    }

    public void setMailing_name(String mailing_name) {
        this.mailing_name = mailing_name;
    }

    public String getMain_url() {
        return main_url;
    }

    public void setMain_url(String main_url) {
        this.main_url = main_url;
    }
    
    
}
