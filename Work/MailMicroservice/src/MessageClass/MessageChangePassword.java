package MessageClass;

import MyInterface.InterfaceMessages;

public class MessageChangePassword implements InterfaceMessages {
    String PasswordKey = "";
    String toEmail = "";
    String fromEmail = "";
    
    public MessageChangePassword(String PasswordKey, String toEmail, String fromEmail) {
        this.PasswordKey = PasswordKey;
        this.toEmail = toEmail;
        this.fromEmail = fromEmail;
    }

    public String getPasswordKey() {
        return PasswordKey;
    }

    public void setPasswordKey(String PasswordKey) {
        this.PasswordKey = PasswordKey;
    }
    
    @Override
    public String getFromEmail() {
        return fromEmail;
    }

    public void setFromEmail(String fromEmail) {
        this.fromEmail = fromEmail;
    }
    
    @Override
    public String getToEmail() {
        return toEmail;
    }

    public void setToEmail(String toEmail) {
        this.toEmail = toEmail;
    }
    
    @Override
    public String getSubject() {
        return "Смена пароля на сайте utp1.ru";
    }

    @Override
    public String getContent() {
        return "Запрос восстановления пароля<br><br>Пришёл запрос на востановления пароля. Если это были не вы, проигнорируйте данное письмо!<br><br>"
                + "Для получения нового пароля введите код подтверждения:<br><br>" + PasswordKey + "<br>";
    }
}
