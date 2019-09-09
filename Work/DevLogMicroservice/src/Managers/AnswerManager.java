package Managers;

import MyEnum.DevLogEnum;
import MyThread.WriteThread;
import RestObject.*;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.rabbitmq.client.AMQP;
import com.rabbitmq.client.Channel;
import com.rabbitmq.client.DefaultConsumer;
import com.rabbitmq.client.Envelope;
import java.io.IOException;

public class AnswerManager extends DefaultConsumer {
    private Channel channel;
    Gson gson;
    WriteThread myThread;
    MessageManager message_manager;
    LogManager log_manager;
    RamManager ram_manager;
    
    String name_class = "AnswerManagerDevLog";
    
    public AnswerManager(Channel channel, GlobalVariables global_variables, LogManager log_manager) {
        super(channel);
        this.channel = channel;
        this.log_manager = log_manager;
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());  
            builder.registerTypeAdapter(DevLog.class, new DevLogConverter());
            builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
            gson = builder.create();
            
            ram_manager = new RamManager(log_manager, message_manager, gson);
            message_manager = new MessageManager(global_variables);
            
            myThread = new WriteThread();
            Thread thread = new Thread(myThread, "My thread");
            thread.start();
        } catch (Exception ex) {
            log_manager.CallError(ex.toString(), "AnswerManagerDevLog. Error on init", name_class);
        } 
    }
    
    @Override
    public void handleDelivery(String string, Envelope envlp, AMQP.BasicProperties bp, byte[] bytes) throws IOException {    
        AMQP.BasicProperties replyProps = new AMQP.BasicProperties.Builder()
        .correlationId(bp.getCorrelationId())
        .build();
        
        String response = new String(bytes, "UTF-8");
        RestObject rest_obj;
        SQLResult res;

        try {   
            rest_obj = gson.fromJson(response, RestObject.class);
            
            log_manager.CallEvent(" receive message: " + rest_obj.toString(DevLogEnum.values()[rest_obj.getCode()].name()), name_class);
            System.out.println("receive message: " + rest_obj.toString(DevLogEnum.values()[rest_obj.getCode()].name()));
            
            res = ram_manager.CheckSessionKey(rest_obj.getSession_key());
            if (res.getId() == -1 || !"Ok".equals(res.getResult())) {
                response = gson.toJson(new SQLResult("auth-failed"), SQLResult.class);
                
            } 
            else {
                rest_obj.getJson_message().addProperty("user_id", res.getId());
                switch (DevLogEnum.values()[rest_obj.getCode()]) {
                    case Write_log:
                        myThread.setReportsString(gson.fromJson(rest_obj.getJson_message().toString(), DevLog.class).toString());
                        response = gson.toJson(new SQLResult("Ok"), SQLResult.class);                    
                        break;
                }
            }
        }
        catch(Exception e) {
            response = log_manager.CallError(e.toString(), "handleDelivery(); Error to send.", name_class);
        }
        finally {
            log_manager.CallEvent("send message: " + response, name_class);
            System.out.println("send message: " + response);
            channel.basicPublish("", bp.getReplyTo(), replyProps, response.getBytes("UTF-8"));
            channel.basicAck(envlp.getDeliveryTag(), false);
        }
    }
}
