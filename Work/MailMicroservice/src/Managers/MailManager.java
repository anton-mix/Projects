package Managers;

import java.util.Properties;
import javax.mail.Message;
import javax.mail.MessagingException;
import javax.mail.PasswordAuthentication;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeMessage;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import SendAndResivClass.*;
import MyInterface.InterfaceMessages;
import MessageClass.MessageAuth;
import MessageClass.MessageAuthConverter;

public class MailManager {   
    String username = "root@utp1.ru"; 
    String password = "WPW6kHujlxOmUG0RmaFsPk2rQksg8F2f";
    String host = "utp1.ru";
    LogManager log_manager;
    Gson gson;

    public MailManager(GlobalVariables global_var) {
        log_manager = new LogManager(global_var, "MailManager");
        this.username = global_var.getUsername();
        this.password = global_var.getPassword();
        this.host = global_var.getHostMail();

        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(MessageAuth.class, new MessageAuthConverter());  
        builder.registerTypeAdapter(Result.class, new ResultConverter());  
        gson = builder.create();      
    }
    
    public <T extends InterfaceMessages>  String SendMineMessage(T item) {
        Properties props = new Properties(); 
        props.put("mail.smtps.user", username);
        props.put("mail.smtsp.password", password);
        props.put("mail.transport.protocol", "smtps");
        props.put("mail.smtps.host", host);
        props.put("mail.smtps.port", "465");
        props.put("mail.smtps.starttls.enable", "true");

        Session session = Session.getDefaultInstance(props, null);

        try {
            Message message = new MimeMessage(session);
            message.setFrom(new InternetAddress(item.getFromEmail()));
            //message.setFrom(new InternetAddress("ipk.stalin@yandex.ru"));
            
            message.setReplyTo(InternetAddress.parse(item.getFromEmail()));
            message.setRecipients(Message.RecipientType.TO, InternetAddress.parse(item.getToEmail()));
            message.setSubject(item.getSubject());
            message.setContent(item.getContent(), "text/html;charset=UTF-8");
            Transport.send(message);
        } catch (MessagingException ex) {
            return log_manager.CallError(ex.toString(), "SendMineMessage(String, String). Exeption In send message");
        }

        return gson.toJson(new Result("Ok"), Result.class);
    }
}