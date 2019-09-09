package Managers;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.rabbitmq.client.AMQP;
import com.rabbitmq.client.Channel;
import com.rabbitmq.client.DefaultConsumer;
import com.rabbitmq.client.Envelope;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.naming.InitialContext;
import javax.naming.NamingException;
import javax.sql.DataSource;

import MessageClass.*;
import MyEnum.*;
import SendAndResivClass.*;

public class AnswerManager extends DefaultConsumer {
    LogManager log_manager;
    GlobalVariables global_variables;
    Channel channel;
    Gson gson;
    
    MailManager mail_manager;

    public AnswerManager(Channel channel, GlobalVariables global_variables) {
        super(channel);
        this.channel = channel;
        this.global_variables = global_variables;
        
        mail_manager = new MailManager(global_variables);
        log_manager = new LogManager(global_variables, "AnswerManagerMail");   

        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(Result.class, new ResultConverter());
        builder.registerTypeAdapter(RessiveClass.class, new RessiveClassConverter());
        builder.registerTypeAdapter(MessageModel.class, new MessageModelConverter());
        builder.registerTypeAdapter(MessageAuth.class, new MessageAuthConverter());
        builder.registerTypeAdapter(MessageChangePassword.class, new MessageChangePasswordConverter());
        builder.registerTypeAdapter(SimpleMessage.class, new SimpleMessageConverter());
        gson = builder.create();
    }

    @Override
    public void handleDelivery(String string, Envelope envlp, AMQP.BasicProperties bp, byte[] bytes) throws IOException {          
        AMQP.BasicProperties replyProps = new AMQP.BasicProperties.Builder()
                .correlationId(bp.getCorrelationId())
                .build();

        String response = new String(bytes, "UTF-8");

        log_manager.setErrorHttp(response);

        try {
            RessiveClass ressiv_obj = gson.fromJson(response, RessiveClass.class);
            System.out.println(" receive message: " + ressiv_obj.toString());
            log_manager.CallEvent(" receive message: " + ressiv_obj.toString());   

            switch(MailMessageEnum.values()[ressiv_obj.getCode()]) {
                case SendMail:
                    response = mail_manager.SendMineMessage(gson.fromJson(ressiv_obj.getJson_message(), MessageAuth.class));
                break;
                case SendChangePassword:
                    response = mail_manager.SendMineMessage(gson.fromJson(ressiv_obj.getJson_message(), MessageChangePassword.class));
                break;
                case SendOtherMessage:
                    response = mail_manager.SendMineMessage(gson.fromJson(ressiv_obj.getJson_message(), SimpleMessage.class));
                break;

            }                    
        } 
        catch (Exception e) {
            response = log_manager.CallError(e.toString(), "call to queue");
        } 
        finally {

            if(bp.getReplyTo() == null) {
                log_manager.CallEvent("send message: " + response + " is drop");
                channel.basicAck(envlp.getDeliveryTag(), false);
            }
            else {
                log_manager.CallEvent("send message: " + response);
                System.out.println("send message: " + response);
                channel.basicPublish("", bp.getReplyTo(), replyProps, response.getBytes("UTF-8"));
                channel.basicAck(envlp.getDeliveryTag(), false);
            }
            synchronized (this) {
                this.notify();
            }
        }
    }
}