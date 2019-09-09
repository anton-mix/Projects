package mailmicroservice;

import Managers.AnswerManager;
import Managers.LogManager;
import SendAndResivClass.GlobalVariables;

import com.rabbitmq.client.ConnectionFactory;
import com.rabbitmq.client.Connection;
import com.rabbitmq.client.Channel;
import com.rabbitmq.client.Consumer;

import java.io.IOException;
import java.util.concurrent.TimeoutException;

public class MailMicroservice {
    
    static LogManager log_manager = null;
    static GlobalVariables global_variables;
    
    public static void main(String[] args) {
        global_variables = new GlobalVariables(args[0], args[1], args[2], Integer.valueOf(args[3]), args[4], args[5], args[6]);
        log_manager = new LogManager(global_variables, "MailMicroservice");
        
        Connection connection = null;
        ConnectionFactory factory = null;
        try {  
            factory = new ConnectionFactory();
            factory.setHost(global_variables.getHost());
            factory.setPort(5672);
            connection = factory.newConnection();
            
            final Channel channel = connection.createChannel();
            channel.queueDeclare(global_variables.getQueueName(), false, false, false, null);
            channel.basicQos(1);
            
            System.out.println(" [x] Host: " + global_variables.getHost() + ", Queue: " + global_variables.getQueueName());
            
            Consumer consumer = new AnswerManager(channel, global_variables);

            channel.basicConsume(global_variables.getQueueName(), false, consumer);
            
            log_manager.CallEvent("Start");
            
            while (true) {
                synchronized(consumer) {
                    try {
      			consumer.wait();
                    } catch (InterruptedException e) {
                        log_manager.CallError(e.toString(), "main. Main wait");	    	
                    }
                }
            }
        } 
        catch (IOException | TimeoutException e) {
            log_manager.CallError(e.toString(), "main. Main program");
        }
        finally {
            if (connection != null)
                try {
                    connection.close();
                } 
                catch (IOException _ignore) {
                    log_manager.CallError(_ignore.toString(), "main. Close connection");
                }
        }
    }
}