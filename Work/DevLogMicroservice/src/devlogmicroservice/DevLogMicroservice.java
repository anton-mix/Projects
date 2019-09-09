
package devlogmicroservice;

import Managers.AnswerManager;
import Managers.LogManager;
import RestObject.GlobalVariables;
import com.rabbitmq.client.*;

public class DevLogMicroservice {
    static LogManager log_manager = null;
    static GlobalVariables global_variables;
    static String name_class = "DevLogMicroservice";
    
    public static void main(String[] argv) throws Exception {
        global_variables = new GlobalVariables(argv[0], argv[1], argv[2], Integer.valueOf(argv[3]));
        log_manager = new LogManager(global_variables, "DevLogMicroservice");
        
        ConnectionFactory factory = new ConnectionFactory();
        factory.setHost(global_variables.getHost());
        factory.setPort(global_variables.getPort());
                
        final Connection connection = factory.newConnection();
        final Channel channel = connection.createChannel();

        channel.queueDeclare(global_variables.getQueue_name(), true, false, false, null);
        
        System.out.println(" [*] Waiting for messages. To exit press CTRL+C");

        channel.basicQos(1);
        final Consumer consumer = new AnswerManager(channel, global_variables, log_manager);
        
        channel.basicConsume(global_variables.getQueue_name(), false, consumer);

        log_manager.CallEvent("Start", name_class);

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
