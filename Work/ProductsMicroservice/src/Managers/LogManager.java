package Managers;

import java.io.IOException;
import java.util.concurrent.TimeoutException;
import com.rabbitmq.client.ConnectionFactory;
import com.rabbitmq.client.Connection;
import com.rabbitmq.client.Channel;
import com.rabbitmq.client.MessageProperties;
import java.util.Date;

public class LogManager {
    String queue_name = "", host = "", name_class;

    public LogManager(String host, String queueName, String name_class) {
        this.queue_name = queueName;
        this.host = host;
        this.name_class = name_class;
    }
    
    public Boolean Call(String Message) {
        try {
            ConnectionFactory factory = new ConnectionFactory();
            factory.setHost(host);
            Connection connection = factory.newConnection();
            Channel channel = connection.createChannel();

            channel.queueDeclare(queue_name, true, false, false, null);

            channel.basicPublish( "", queue_name, MessageProperties.PERSISTENT_TEXT_PLAIN, Message.getBytes());

            channel.close();
            connection.close();
            return true;
        }
        catch(IOException | TimeoutException ex) {
            return false;   
        }
    }
    
    public String CallEvent (String event) {
        if(Call("Type: event, Time: " + new Date().toString() + ", Name: " + name_class + ", Event: " + event)) {
            return "[{\"Result\":\"Ok\" }]";
        } else {
            System.out.println("Error log");
            return "[{\"Result\":\"ErrorLog\" }]";
        }
    }   
    
    public String CallError (String TextError, String event) {
        if(Call("Type: error, Time: " + new Date().toString() + ", Name: " + name_class + ", Event: Error in " + event + ", Error: " + TextError )) {
            System.out.println("Error");
            return "[{\"Result\":\"Error\" }]";
        } else {
            System.out.println("Error log");
            return "[{\"Result\":\"ErrorLog\" }]";
        }
    }
}
