/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Microservice;

import Managers.AnswerManager;
import Managers.ProductsManager;
import Managers.LogManager;
import SQLClass.GlobalVariables;
import com.rabbitmq.client.Channel;
import com.rabbitmq.client.Connection;
import com.rabbitmq.client.ConnectionFactory;
import com.rabbitmq.client.Consumer;
import java.io.IOException;
import java.util.concurrent.TimeoutException;

public class ProductsMicroservice {
    static LogManager log_manager = null;
    static ProductsManager products_manager;
    static GlobalVariables global_variables;
    
    public static void main(String[] args) {
        global_variables = new GlobalVariables(args[0], args[1], args[2], args[3], args[4], args[5], args[6]);
        log_manager = new LogManager(global_variables.getHost(), global_variables.getQueueNameLog(), "ProductsMicroservice");
        
        ConnectionFactory factory = new ConnectionFactory();
        factory.setHost(global_variables.getHost());

        Connection connection = null;
        try {
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
            log_manager.CallError(e.toString(), "main. Main programn");
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
