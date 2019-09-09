package Managers;

import SendAndResivClass.*;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.io.IOException;
import com.rabbitmq.client.ConnectionFactory;
import com.rabbitmq.client.Connection;
import com.rabbitmq.client.Channel;
import com.rabbitmq.client.MessageProperties;
import java.util.Date;
import java.util.concurrent.TimeoutException;

public class LogManager {
    String name_class = "";
    String ErrorHttp = "";
    Gson gson;
    GlobalVariables global_variables;

    public LogManager(GlobalVariables global_variables, String name_class) {
        this.global_variables = global_variables;
        this.name_class = name_class;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(Result.class, new ResultConverter());
        gson = builder.create();
    }
    
    public Boolean Call(String Message) {
        try {
            ConnectionFactory factory = new ConnectionFactory();
            factory.setHost(global_variables.getHost());
            factory.setPort(global_variables.getPort());
            Connection connection = factory.newConnection();
            Channel channel = connection.createChannel();

            channel.queueDeclare(global_variables.getQueueNameLog(), true, false, false, null);

            channel.basicPublish( "", global_variables.getQueueNameLog(), MessageProperties.PERSISTENT_TEXT_PLAIN, Message.getBytes());

            channel.close();
            connection.close();
            return true;
        }
        catch(IOException | TimeoutException ex) {
            return false;   
        }
    }
    
    public String CallEvent (String event) {
        if(Call("Type: event, Time: " + new Date().toString() + ", Name: " + name_class + ", Event:" + event)) {
            return gson.toJson(new Result("Ok"));
        } else {
            System.out.println("Error log");
            return gson.toJson(new Result("ErrorLog"));
        }
    }   
    
    public String CallError (String TextError, String event) {
        Call("Type: error, Time: " + new Date().toString() + ", Name: " + name_class + ", Event: Error in " + event + ", Error: " + ErrorHttp );
        if(Call("Type: error, Time: " + new Date().toString() + ", Name: " + name_class + ", Event: Error in " + event + ", Error: " + TextError )) {
            System.out.println("Error");
            return gson.toJson(new Result("Error"));
        } else {
            System.out.println("Error log");
            return gson.toJson(new Result("ErrorLog"));
        }
    }

    public String getErrorHttp() {
        return ErrorHttp;
    }

    public void setErrorHttp(String ErrorHttp) {
        this.ErrorHttp = ErrorHttp;
    }
}
