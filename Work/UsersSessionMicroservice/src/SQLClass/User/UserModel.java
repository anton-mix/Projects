package SQLClass.User;

import MyInterface.ModelInterface;
import QueueClass.UserAuthModel;
import java.util.ArrayList;

public class UserModel implements ModelInterface {
    int id = -1;
    String username = "";
    String password = "";
    String surname = "";
    String name = "";
    String patronymic = "";
    String avatar_url = "";
    String email = "";
    String phone = "";

    public UserModel(int id, String username, String password, String surname, 
            String name, String gpatronymic, String avatar_url, String email, String phone) {
        this.id = id;
        this.username = username;
        this.password = password;
        this.surname = surname;
        this.name = name;
        this.patronymic = gpatronymic;
        this.avatar_url = avatar_url;
        this.email = email;
        this.phone = phone;
    }
    
    public UserModel(String username, String password, String surname, 
            String name, String gpatronymic, String avatar_url, String email, String phone) {
        this.username = username;
        this.password = password;
        this.surname = surname;
        this.name = name;
        this.patronymic = gpatronymic;
        this.avatar_url = avatar_url;
        this.email = email;
        this.phone = phone;
    }
        
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), username, password, surname, name, patronymic, avatar_url, email, phone};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            "'" + username.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + password.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + surname.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + patronymic.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + avatar_url.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + email.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + phone.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "username", "password", "surname", "name", "gpatronymic", "avatar_url", "email", "phone"}; 
    }
    
    @Override
    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
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

    public String getSurname() {
        return surname;
    }

    public void setSurname(String surname) {
        this.surname = surname;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getGpatronymic() {
        return patronymic;
    }

    public void setGpatronymic(String gpatronymic) {
        this.patronymic = gpatronymic;
    }

    public String getAvatar_url() {
        return avatar_url;
    }

    public void setAvatar_url(String avatar_url) {
        this.avatar_url = avatar_url;
    }

    public String getPatronymic() {
        return patronymic;
    }

    public void setPatronymic(String patronymic) {
        this.patronymic = patronymic;
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
    
    public void addFilds(UserAuthModel item) {
        this.username = item.getUsername();
        this.password = item.getPassword();
        this.name = item.getName();
        this.email = item.getEmail(); 
        this.phone = item.getPhone();
    }
}

