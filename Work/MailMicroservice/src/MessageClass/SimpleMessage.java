package MessageClass;

import MyInterface.InterfaceMessages;

public class SimpleMessage implements InterfaceMessages {
    String text = "";
    String toEmail = "";
    String fromEmail = "";
    
    public SimpleMessage(String text, String toEmail, String fromEmail) {
        this.text = text;
        this.toEmail = toEmail;
        this.fromEmail = fromEmail;
    }

    public String getText() {
        return text;
    }

    public void setText(String text) {
        this.text = text;
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
        return "Пробное сообщение";
    }

    @Override
    public String getContent() {
        return text;
    }
}

