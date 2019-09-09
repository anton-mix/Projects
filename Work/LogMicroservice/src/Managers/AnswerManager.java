package Managers;

import MyThreads.WriteThread;
import RestObject.GlobalVariables;
import RestObject.RestObject;
import RestObject.RestObjectConverter;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.rabbitmq.client.AMQP;
import com.rabbitmq.client.Channel;
import com.rabbitmq.client.DefaultConsumer;
import com.rabbitmq.client.Envelope;
import java.io.IOException;
import java.util.concurrent.TimeoutException;
import java.util.logging.Level;
import java.util.logging.Logger;

public class AnswerManager  extends DefaultConsumer {
    private Channel channel;
    Gson gson;
    WriteThread myThread;
    MessageManager message_manager;
    String name_class = "AnswerManagerLog";
    
    public AnswerManager(Channel channel, GlobalVariables global_variables) {
        super(channel);
        this.channel = channel;
            
         try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());  
            builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
            gson = builder.create();

            message_manager = new MessageManager(global_variables);
            myThread = new WriteThread(message_manager, gson);
            Thread thread = new Thread(myThread, "My thread");
            thread.start();
         } catch (IOException ex) {
             Logger.getLogger(AnswerManager.class.getName()).log(Level.SEVERE, null, ex);
         } catch (TimeoutException ex) {
             Logger.getLogger(AnswerManager.class.getName()).log(Level.SEVERE, null, ex);
         }  catch (Exception ex) {
             Logger.getLogger(AnswerManager.class.getName()).log(Level.SEVERE, null, ex);
         }  
    }
    
    @Override
    public void handleDelivery(String string, Envelope envlp, AMQP.BasicProperties bp, byte[] bytes) throws IOException {    
        String response = new String(bytes, "UTF-8");
        
        try {
            System.out.println(" [x] Received '" + response + "'");
            String buff = response.split(" ")[1];

            if ("event,".equals(buff)) {
                myThread.setEventString(response);
            }
            if ("error,".equals(buff)) {
                myThread.setErrorString(response);
            }
        }  catch (Exception ex) {
            System.out.println(ex.toString());
        }  finally {
            channel.basicAck(envlp.getDeliveryTag(), false);
        }
    }
}

