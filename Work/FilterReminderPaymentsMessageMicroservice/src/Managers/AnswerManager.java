package Managers;

import MyEnum.FilterReminderPaymentsMessageEnum;
import MyInterface.InterfaceConverter;
import RestObject.RestObject;
import RestObject.RestObjectConverter;
import SQLClass.GlobalVariables;
import SQLClass.SQLResult;
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

public class AnswerManager  extends DefaultConsumer implements InterfaceConverter {
    private Channel channel;
    Connection connection_with_db = null;
    LogManager log_manager;
    MessageManager message_manager;
    MemCachedManager memCached_manager;
    Gson gson;
    
    String connection_string;   
    String name_class = "FilterSingleMessageAnswerManager";
       
    public AnswerManager(Channel channel, GlobalVariables global_variables, LogManager log_manager) {
        super(channel);
        this.channel = channel;
        this.log_manager = log_manager;
   
        try {
            Class.forName("com.mysql.jdbc.Driver").newInstance();
            connection_string = "jdbc:mysql://" + global_variables.getDb_ip_sql() + "/" + global_variables.getDb_name() + "?useSSL=false&user=" + global_variables.getDb_user_name() + "&password=" + global_variables.getDb_password() + "&useUnicode=true&characterEncoding=UTF-8&serverTimezone=UTC";             
            
            message_manager = new MessageManager(global_variables);
            memCached_manager = new MemCachedManager(log_manager, global_variables);
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter()); 
            gson = builder.create();
     
        } catch (Exception ex) {
            log_manager.CallError(ex.toString(), "AnswerManager. Error on init", name_class);
        } 
    }

    @Override
    public void handleDelivery(String string, Envelope envlp, AMQP.BasicProperties bp, byte[] bytes) throws IOException {     
        AMQP.BasicProperties replyProps = new AMQP.BasicProperties.Builder()
            .correlationId(bp.getCorrelationId())
            .build();
        
        String response = new String(bytes, "UTF-8");
        SQLResult res;
            
        log_manager.setErrorHttp(response);
        
        try {
            RestObject rest_obj = gson.fromJson(response, RestObject.class);
            
            System.out.println("receive message: " + rest_obj.toString(FilterReminderPaymentsMessageEnum.values()[rest_obj.getCode()].name()));
            log_manager.CallEvent(" receive message: " + rest_obj.toString(FilterReminderPaymentsMessageEnum.values()[rest_obj.getCode()].name()), name_class);
            
            connection_with_db = DriverManager.getConnection(connection_string);
            
            log_manager.setErrorHttp(FilterReminderPaymentsMessageEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
            
            switch (FilterReminderPaymentsMessageEnum.values()[rest_obj.getCode()]) {
                case tableRefresh:
                    break;
                case addItems:
                    break;
            }
        }
        catch(Exception e) {
            log_manager.CallError(e.toString(), "handleDelivery(); Error to send.", name_class);
        }
        finally {
            
            try {
                connection_with_db.close();
            } catch (SQLException ex) {
                log_manager.CallEvent("handleDelivery(); Connection Close Error", name_class);
            }
            
            if(bp.getReplyTo() == null) {
                log_manager.CallEvent("send message: " + response + " is drop", name_class);
                System.out.println("send message: " + response + " is drop");
                channel.basicAck(envlp.getDeliveryTag(), false);
            }
            else {
                log_manager.CallEvent("send message: " + response, name_class);
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
