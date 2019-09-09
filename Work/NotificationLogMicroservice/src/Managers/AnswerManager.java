package Managers;

import MyEnum.RestNotificationLogEnum;
import SQLClass.GlobalVariables;
import MyInterface.InterfaceConverter;

import RestObject.RestObject;
import RestObject.RestObjectConverter;
import SQLClass.NotificationLog.NotificationLogModel;
import SQLClass.NotificationLog.NotificationLogModelConverter;
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
import java.sql.Timestamp;

public class AnswerManager extends DefaultConsumer implements InterfaceConverter {
    private Channel channel;
    Connection connection_with_db = null;
    GlobalVariables global_variables;
    Gson gson;
    
    String connection_string;   
    String name_class = "NotificationLogAnswerManager";
    
    MessageManager message_manager;
    NotificationLogManager notification_Log_Manager;
    LogManager log_manager;
    RamManager ram_manager;
    
    public AnswerManager(Channel channel, GlobalVariables global_variables, LogManager log_manager) {
        super(channel);
        this.channel = channel;
        this.global_variables = global_variables;
        this.log_manager = log_manager;
   
        try {
           
            Class.forName("com.mysql.jdbc.Driver").newInstance();
            connection_string = "jdbc:mysql://" + global_variables.getDb_ip_sql()+ "/" + global_variables.getDb_name()+ "?useSSL=false&user=" + global_variables.getDb_user_name() + "&password=" + global_variables.getDb_password() + "&useUnicode=true&characterEncoding=UTF-8&serverTimezone=UTC";             
            
            message_manager = new MessageManager(global_variables);
            notification_Log_Manager = new NotificationLogManager(global_variables, null, log_manager);
                    
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());  
            builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());  
            builder.registerTypeAdapter(NotificationLogModel.class, new NotificationLogModelConverter());
            gson = builder.create();    
            ram_manager = new RamManager(log_manager, message_manager, gson);
        } 
        catch (Exception ex) {
            log_manager.CallError(ex.toString(), "AnswerManagerMailing. Error on init", name_class);
        } 
    }

    @Override
    public void handleDelivery(String string, Envelope envlp, AMQP.BasicProperties bp, byte[] bytes) throws IOException {     
        AMQP.BasicProperties replyProps = new AMQP.BasicProperties.Builder()
                .correlationId(bp.getCorrelationId())
                .build();
        ResData response = new ResData();
        response.setRes(new String(bytes, "UTF-8"));
        SQLResult res;

        try {
            RestObject rest_obj = gson.fromJson(response.getRes(), RestObject.class);
            
            log_manager.CallEvent(" receive message: " + rest_obj.toString(RestNotificationLogEnum.values()[rest_obj.getCode()].name()), name_class);
            System.out.println("receive message: " + rest_obj.toString(RestNotificationLogEnum.values()[rest_obj.getCode()].name()));
            
            connection_with_db = DriverManager.getConnection(connection_string);
            
            notification_Log_Manager.setConnection_with_db(connection_with_db);
            
            log_manager.setErrorHttp(RestNotificationLogEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            notification_Log_Manager.setErrorHttp(RestNotificationLogEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            
            res = ram_manager.CheckSessionKey(rest_obj.getSession_key());
            if (res.getId() == -1 || !"Ok".equals(res.getResult())) {
                response = new ResData(false, "{\"result\":\"auth-failed\"}");
            } else {
                rest_obj.getJson_message().addProperty("user_id", res.getId());
                boolean need_commit = false;
                switch (RestNotificationLogEnum.values()[rest_obj.getCode()]) {
                    case Notification_log_select:
                        response = notification_Log_Manager.getNotificationLog(gson.fromJson(rest_obj.getJson_message(), NotificationLogModel.class));
                        break;
                    case Notification_log_select_count:
                        response = notification_Log_Manager.getCountNotificationLog(gson.fromJson(rest_obj.getJson_message(), NotificationLogModel.class));
                        break;
                    case Notification_log_select_with_count: 
                        response = notification_Log_Manager.Notification_log_select_with_count(gson.fromJson(rest_obj.getJson_message(), NotificationLogModel.class));
                        break;
                    case Notification_log_add:
                        response = notification_Log_Manager.AddNotificationLog(gson.fromJson(rest_obj.getJson_message(), NotificationLogModel.class));
                        need_commit = true;
                        break;
                    case Notification_log_update:
                        response = notification_Log_Manager.UpdateNotificationLog(gson.fromJson(rest_obj.getJson_message(), NotificationLogModel.class));
                        need_commit = true;
                        break;
                    case Notification_log_delete:
                        response = notification_Log_Manager.DeleteNotificationLog(gson.fromJson(rest_obj.getJson_message(), NotificationLogModel.class));
                        need_commit = true;
                        break;
                }
                if (need_commit) {
                    if (response.getIs_success()) {
                        connection_with_db.commit();
                    } else {
                        connection_with_db.rollback();
                    }
                }
            }
        }
        catch(Exception e) {
            response = log_manager.CallError(e.toString(), "handleDelivery(); Error to send.", name_class);
        }
        finally {
            try {
                connection_with_db.close();
            } catch (SQLException ex) {
                log_manager.CallEvent("handleDelivery(); Connection Close Error", name_class);
            }
            
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
