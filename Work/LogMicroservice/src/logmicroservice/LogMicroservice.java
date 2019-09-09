package logmicroservice;

import Managers.AnswerManager;
import MyThreads.WriteThread;
import RestObject.GlobalVariables;
import com.rabbitmq.client.*;

public class LogMicroservice {

    GlobalVariables global_variables;
    static WriteThread myThread;
    
    public static void main(String[] argv) throws Exception {
        GlobalVariables global_variables = new GlobalVariables(argv[0], argv[1], Integer.valueOf(argv[2]));
        
        ConnectionFactory factory = new ConnectionFactory();
        factory.setHost(global_variables.getHost());
        factory.setPort(global_variables.getPort());
                
        final Connection connection = factory.newConnection();
        final Channel channel = connection.createChannel();

        channel.queueDeclare(global_variables.getQueue_name(), true, false, false, null);
        
        System.out.println(" [*] Waiting for messages. To exit press CTRL+C");

        channel.basicQos(1);

        final Consumer consumer = new AnswerManager(channel, global_variables);
        
        channel.basicConsume(global_variables.getQueue_name(), false, consumer);

        while (true) {
            synchronized (consumer) {
                try {
                    consumer.wait();
                } catch (InterruptedException e) {

                }
            }
        }
    }
}
