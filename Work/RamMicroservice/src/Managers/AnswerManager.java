package Managers;

import SQLClass.GlobalVariables;
import MyInterface.InterfaceConverter;
import QueueClass.ArrUserPasswordKeyValue;
import QueueClass.ArrUserPasswordKeyValueConverter;
import QueueClass.ResQueueItem;
import QueueClass.UserAuthModel;
import QueueClass.UserAuthModelConverter;
import RestObject.RestObject;
import RestObject.RestObjectConverter;

import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.rabbitmq.client.AMQP;
import com.rabbitmq.client.Channel;
import com.rabbitmq.client.DefaultConsumer;
import com.rabbitmq.client.Envelope;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Date;
import java.util.Map;

public class AnswerManager extends DefaultConsumer implements InterfaceConverter {
    private Channel channel;
    GlobalVariables global_variables;
    LogManager log_manager;
    Gson gson;
   
    String name_class = "RamAnswerManager";
    
    SessionKeyQueueManager SessionKeyQueue_manager;
    UserAuthMessageQueueManager UserAuthMessageQueue_manager;
    ChangePassMessageQueueManager ChangePassMessageQueue_manager;
       
    public AnswerManager(Channel channel, GlobalVariables global_variables, LogManager log_manager) {
        super(channel);
        this.channel = channel;
        this.global_variables = global_variables;
        
        this.log_manager = log_manager;
   
        try { 
            
            SessionKeyQueue_manager = new SessionKeyQueueManager(log_manager);
            UserAuthMessageQueue_manager = new UserAuthMessageQueueManager(log_manager);
            ChangePassMessageQueue_manager = new ChangePassMessageQueueManager(log_manager);
    
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            builder.registerTypeAdapter(UserAuthModel.class, new UserAuthModelConverter());
            builder.registerTypeAdapter(ArrUserPasswordKeyValue.class, new ArrUserPasswordKeyValueConverter());
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
        
        ResData response = new ResData();
        response.setRes(new String(bytes, "UTF-8"));
            
        log_manager.setErrorHttp(response.getRes());
        
        try {
            RestObject rest_obj = gson.fromJson(response.getRes(), RestObject.class);
            
            System.out.println("receive message: " + rest_obj.toString());
            log_manager.CallEvent(" receive message: " + rest_obj.toString(), name_class);

            log_manager.setErrorHttp(rest_obj.getCode().name()  + " " + rest_obj.getJson_message().toString());
           
            switch (rest_obj.getCode()) {
                case Get_session_key:
                    response = SessionKeyQueue_manager.GetSessionKeyItem(gson.fromJson(rest_obj.getJson_message(), ResQueueItem.class));
                    break;
                case Add_session_key:
                    response = SessionKeyQueue_manager.addSessionKeyItem(gson.fromJson(rest_obj.getJson_message(), ResQueueItem.class));
                    break;
                case Get_session_keys:
                    response = SessionKeyQueue_manager.GetSessionKeyItems();
                    break;
                case Check_session_key:
                    response = SessionKeyQueue_manager.CheckUserId(rest_obj.getSession_key());
                    break;
                case Add_auth_message_key:
                    response = UserAuthMessageQueue_manager.AddUserAuthMessageItem(gson.fromJson(rest_obj.getJson_message(), UserAuthModel.class));
                    break;
                case Get_auth_message_keys:
                    response = UserAuthMessageQueue_manager.GetUserAuthMessageItems();
                    break;
                case Check_auth_message_key:
                    response = UserAuthMessageQueue_manager.CheckAuthKey(rest_obj.getSession_key());
                    break;
                case Delete_auth_message_key:
                    response = UserAuthMessageQueue_manager.DeleteUserAuthMessageItem(rest_obj.getSession_key());
                    break;
                case Delete_user_auth_message_item_for_email:
                    response = UserAuthMessageQueue_manager.DeleteUserAuthMessageItemForEmail(gson.fromJson(rest_obj.getJson_message(), UserAuthModel.class));
                    break;
                case Check_and_delete_user_auth_message:
                    response = UserAuthMessageQueue_manager.CheckAndDeleteUserAuthMessage(rest_obj.getSession_key());
                    break;
                case Check_auth_message_key_email:
                    response = UserAuthMessageQueue_manager.CheckAuthMessageKeyEmail(rest_obj.getSession_key());
                    break;
                case Add_change_pass_message_key:
                    response = ChangePassMessageQueue_manager.AddChangePassMessageItem(gson.fromJson(rest_obj.getJson_message(), ArrUserPasswordKeyValue.class));
                    break;
                case Get_change_pass_message_keys:
                    response = ChangePassMessageQueue_manager.GetChangePassMessageItems();
                    break;
                case Check_change_pass_message_key:
                    response = ChangePassMessageQueue_manager.CheckVhangePassKey(rest_obj.getSession_key());
                    break;
                case Check_change_pass_message_key_login:
                    response = ChangePassMessageQueue_manager.CheckChangePassMessageKeyLogin(rest_obj.getSession_key());
                    break;
                case Delete_change_pass_message_key:
                    response = ChangePassMessageQueue_manager.DeleteChangePassMessageItem(rest_obj.getSession_key());
                    break;
                case Delete_change_pass_message_item_from_login:
                    response = ChangePassMessageQueue_manager.DeleteChangePassMessageItemFromLogin(gson.fromJson(rest_obj.getJson_message(), ArrUserPasswordKeyValue.class));
                    break;
                case Check_and_delete_change_pass_message_key:
                    response = ChangePassMessageQueue_manager.CheckAndDeleteChangePassMessageKey(rest_obj.getSession_key());
                    break;
                case Refres_queues:
                    Date dateNow = new Date();
                    response = SessionKeyQueue_manager.RefresQueue(dateNow, 3600000);
                    if(!response.getIs_success())
                        break;
                    response = UserAuthMessageQueue_manager.RefresQueue(dateNow, 3600000);
                    if (!response.getIs_success()) {
                        break;
                    }
                    response = ChangePassMessageQueue_manager.RefresQueue(dateNow, 3600000);
                    if (!response.getIs_success()) {
                        break;
                    }
                    break;
            }
        }
        catch(Exception e) {
            log_manager.CallError(e.toString(), "handleDelivery(); Error to send.", name_class);
        }
        finally {    
            if(bp.getReplyTo() == null) {
                log_manager.CallEvent("send message: " + response.getRes() + " is drop", name_class);
                System.out.println("send message: " + response.getRes() + " is drop");
                channel.basicAck(envlp.getDeliveryTag(), false);
            }
            else {
                log_manager.CallEvent("send message: " + response.getRes(), name_class);
                System.out.println("send message: " + response.getRes());
                channel.basicPublish("", bp.getReplyTo(), replyProps, response.getRes().getBytes("UTF-8"));
                channel.basicAck(envlp.getDeliveryTag(), false);
            }
            
            synchronized (this) {
                this.notify();
            }
        }
    }
}
