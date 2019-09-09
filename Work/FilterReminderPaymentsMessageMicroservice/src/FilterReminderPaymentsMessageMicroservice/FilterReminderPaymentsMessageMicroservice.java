package FilterReminderPaymentsMessageMicroservice;

import Managers.LogManager;
import Managers.MemCachedManager;
import SQLClass.GlobalVariables;
import com.rabbitmq.client.Connection;
import com.rabbitmq.client.ConnectionFactory;
import com.rabbitmq.client.Consumer;

public class FilterReminderPaymentsMessageMicroservice {

    static LogManager log_manager = null;
    static GlobalVariables global_variables;
    
    public static void main(String[] args) {
        global_variables = new GlobalVariables(args[0], args[1], args[2], Integer.valueOf(args[3]), args[4], args[5], args[6], args[7], args[8], Integer.valueOf(args[9]), args[10], args[11]);
        System.out.println();
        log_manager = new LogManager(global_variables, "FilterSingleMessageMicroservice");
        Connection connection = null;
        ConnectionFactory factory = null;
        Consumer consumer = null;
        String name_class = "FilterSingleMessageMicroservice";
        
        MemCachedManager memCached_manager = new MemCachedManager(log_manager, global_variables);
        memCached_manager.SetItem("123", "321");
        System.out.println(memCached_manager.GetItem("123"));
        /*
        try {  
            factory = new ConnectionFactory();
            factory.setHost(global_variables.getHost());
            factory.setPort(global_variables.getPort());
            connection = factory.newConnection();
            
            final Channel channel = connection.createChannel();
            channel.queueDeclare(global_variables.getQueue_name() + "_2", false, false, false, null);
            channel.basicQos(1);

            System.out.println(" [x] Host: " + global_variables.getHost() + ", Queue: " + global_variables.getQueue_name());

            consumer = new AnswerManager(channel, global_variables, log_manager);

            channel.basicConsume(global_variables.getQueue_name(), false, consumer);
            log_manager.CallEvent("Start", name_class);

            while (true) {
                synchronized (consumer) {
                    try {
                        consumer.wait();
                    } catch (InterruptedException e) {
                        log_manager.CallError(e.toString(), "main. Main wait", name_class);
                    }
                }
            }
        } 
        catch (IOException | TimeoutException e) {
            log_manager.CallError(e.toString(), "main. Main program", name_class);
        }
        finally {
            if (connection != null) {
                try {
                    connection.close();
                } catch (Exception ex) {
                    log_manager.CallError(ex.toString(), "main. Close connection", name_class);
                }
            }
        }*/
    }
    
}