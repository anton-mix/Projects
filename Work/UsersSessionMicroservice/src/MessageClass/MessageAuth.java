package MessageClass;

import MyInterface.InterfaceMessages;

public class MessageAuth implements InterfaceMessages{
    String name = "";
    String session_id = "";
    String toEmail = "";
    String fromEmail = "";
    
    public MessageAuth(String name, String session_id, String toEmail, String fromEmail) {
        this.name = name;
        this.session_id = session_id;
        this.toEmail = toEmail;
        this.fromEmail = fromEmail;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getSession_id() {
        return session_id;
    }

    public void setSession_id(String session_id) {
        this.session_id = session_id;
    }

    @Override
    public String getToEmail() {
        return toEmail;
    }

    public void setToEmail(String toEmail) {
        this.toEmail = toEmail;
    }
    
    @Override
    public String getFromEmail() {
        return fromEmail;
    }

    public void setFromEmail(String fromEmail) {
        this.fromEmail = fromEmail;
    }

    @Override
    public String getSubject() {
        return "Подтверждение регистрации на сайте utp1.ru";
    }

    
    
    @Override
    public String getContent() {
        return "Поздравляем с регистрацией на utp1.ru<br><br>Чтобы подтвердить регистрацию, пожалуйста, откройте ссылку:<br>"
                + "<br><a href=\"http://localhost:4200/activation/" + session_id + "\" >https://utp1.ru/activation/" + session_id + "</a>";
    }
}
