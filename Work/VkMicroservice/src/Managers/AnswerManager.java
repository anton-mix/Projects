package Managers;

import MyEnum.RestVkEnum;
import SQLClass.GlobalVariables;
import MyInterface.InterfaceConverter;
import RestObject.RestObject;
import RestObject.RestObjectConverter;
import SQLClass.ResData;


import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import SQLClass.VKAuthClass.VKAuthClass;
import SQLClass.VKAuthClass.VKAuthClassConverter;
import SQLClass.VKSendingClass.VKSendingClass;
import SQLClass.VKSendingClass.VKSendingClassConverter;
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

public class AnswerManager extends DefaultConsumer implements InterfaceConverter {
    private Channel channel;
    GlobalVariables global_variables;
    LogManager log_manager;
    MessageManager message_manager;
    Gson gson;
    
    String connection_string;   
    String name_class = "VKAnswerManager";
    RamManager ram_manager;
    VKManager vk_manager;
       
    public AnswerManager(Channel channel, GlobalVariables global_variables, LogManager log_manager) {
        super(channel);
        this.channel = channel;
        this.global_variables = global_variables;
        
        this.log_manager = log_manager;
   
        try {
            message_manager = new MessageManager(global_variables);
            vk_manager = new VKManager(global_variables, log_manager);
      
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            builder.registerTypeAdapter(VKSendingClass.class, new VKSendingClassConverter()); 
            builder.registerTypeAdapter(VKAuthClass.class, new VKAuthClassConverter());
	    builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
            gson = builder.create();
            
            ram_manager = new RamManager(log_manager, message_manager, gson);
     
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
        SQLResult res;
            
        log_manager.setErrorHttp(response.getRes());
        
        try {
            RestObject rest_obj = gson.fromJson(response.getRes(), RestObject.class);
            System.out.println("receive message: " + rest_obj.toString(RestVkEnum.values()[rest_obj.getCode()].name()));
            log_manager.CallEvent(" receive message: " + rest_obj.toString(RestVkEnum.values()[rest_obj.getCode()].name()), name_class);
            
            log_manager.setErrorHttp(RestVkEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
            vk_manager.setErrorHttp(RestVkEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
            
            //response = message_manager.call(gson.toJson(new RestObject(new JsonObject(), rest_obj.getSession_key(), RestVkEnum.user_check)), "user_session_queue");
            res = ram_manager.CheckSessionKey(rest_obj.getSession_key());
            if (res.getId() == -1 || !"Ok".equals(res.getResult())) {
                response = new ResData(false, "{\"result\":\"auth-failed\"}");
            }
            else {
                rest_obj.getJson_message().addProperty("user_id", res.getId());
                switch (RestVkEnum.values()[rest_obj.getCode()]) {
                    case get_user_groups: 
                        response = vk_manager.getUserGroups(gson.fromJson(rest_obj.getJson_message(), VKAuthClass.class));
                        break;
                        /*
                    case send_messages_by_user: 
                        response = vk_manager.sendMessagesByUser(gson.fromJson(rest_obj.getJson_message(), VKSendingClass.class));
                        break;
                        **/
                    case send_messages_by_group: 
                        System.out.println("aaaaa");
			response = vk_manager.sendMessagesByGroup(gson.fromJson(rest_obj.getJson_message(), VKSendingClass.class));
                        break;
                   /* case get_group_members: 
                        response = vk_manager.getGroupMembers(gson.fromJson(rest_obj.getJson_message(), VKAuthClass.class));
                        break;  **/
                        /*
                    case join_user_to_group:
                        response = vk_manager.joinUserToGroup(gson.fromJson(rest_obj.getJson_message(), VKAuthClass.class));
                        break;
                    case allow_messages_from_group:
                        response = vk_manager.allowMessagesFromGroup(gson.fromJson(rest_obj.getJson_message(), VKAuthClass.class));
                        break;
                        **/
		    default:
			response = new ResData(false,gson.toJson( new SQLResult("invalid arg switch"), SQLResult.class));
			break;
                }
            }
        }
        catch(Exception e) {
            response = log_manager.CallError(e.toString(), "handleDelivery(); Error to send.", name_class);
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
