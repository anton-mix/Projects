package QueueClass;

import MyInterface.ModelInterface;

public class UserAuthModel implements ModelInterface {
    String name = "";
    String username = "";
    String password = "";
    String session_id = "";
    String email = "";
    String phone  = "";

    public UserAuthModel(String name, String username, String password, String session_id, String email, String phone) {
        this.name = name;
        this.username = username;
        this.password = password;
        this.session_id = session_id;
        this.email = email;
        this.phone = phone;
    }
    
    public UserAuthModel(String email) {
        this.name = "";
        this.username = "";
        this.password = "";
        this.session_id = "";
        this.email = email;
        this.phone = "";
    }
    
    @Override
    public String[] getArray() {
        return new String[] {name, username, password, session_id, email, phone};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            "'" + name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + username.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + password.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + session_id.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + email.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + phone.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"name", "username", "password", "session_id", "email", "phone"}; 
    }
    
    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getUsername() {
        return username;
    }

    public void setUsername(String username) {
        this.username = username;
    }

    public String getPassword() {
        return password;
    }

    public void setPassword(String password) {
        this.password = password;
    }

    public String getSession_id() {
        return session_id;
    }

    public void setSession_id(String session_id) {
        this.session_id = session_id;
    }

    public String getEmail() {
        return email;
    }

    public void setEmail(String email) {
        this.email = email;
    }

    public String getPhone() {
        return phone;
    }

    public void setPhone(String phone) {
        this.phone = phone;
    }
    
    public void AddFilds(UserAuthModel item) {
        this.email = ("".equals(item.getEmail()))?this.email:item.getEmail();
        this.name = ("".equals(item.getName()))?this.name:item.getName();
        this.password = ("".equals(item.getPassword()))?this.password:item.getPassword();
        this.phone = ("".equals(item.getPhone()))?this.phone:item.getPhone();
        this.username = ("".equals(item.getUsername()))?this.username:item.getUsername();
    } 
    
    @Override
    public String toString() {
        return " Name: " + name + " username: " + username + " password: " + password + " session_id: " + session_id + " email: " + email + " phone: " + phone;
    }
}
