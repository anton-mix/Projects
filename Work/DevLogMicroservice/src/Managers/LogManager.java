package Managers;

import RestObject.GlobalVariables;
import SQLClass.SQLResultConverter;
import SQLClass.SQLResult;
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
    String ErrorHttp = "";
    String name_microservice = "";
    Gson gson;
    GlobalVariables global_variables;
    
    public LogManager(GlobalVariables global_variables, String name_microservice) {
        this.global_variables = global_variables;
        this.name_microservice = name_microservice;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        gson = builder.create();
    }
    
    public Boolean Call(String Message) {
        try {
            ConnectionFactory factory = new ConnectionFactory();
            factory.setHost(global_variables.getHost());
            factory.setPort(global_variables.getPort());
            Connection connection = factory.newConnection();
            Channel channel = connection.createChannel();

            channel.queueDeclare(global_variables.getQueue_name_log(), true, false, false, null);

            channel.basicPublish( "", global_variables.getQueue_name_log(), MessageProperties.PERSISTENT_TEXT_PLAIN, Message.getBytes());

            channel.close();
            connection.close();
            factory.clone();
            return true;
        }
        catch(IOException | TimeoutException ex) {
            return false;   
        }
    }
    
    public String CallEvent (String event, String name_class) {
        if(Call("Type: event, Time: " + new Date().toString() + ", NameMicroservice: " + name_microservice + ", NameClass: " + name_class + ", Event:" + event)) {
            return gson.toJson(new SQLResult("Ok"));
        } else {
            System.out.println("Error log");
            return gson.toJson(new SQLResult("ErrorLog"));
        }
    }

    public String CallError (String TextError, String event, String name_class) {
        Call("Type: error, Time: " + new Date().toString() + ", NameMicroservice: " + name_microservice + ", NameClass: " + name_class + ", Event: Error in " + event + ", Error: " + ErrorHttp );
        if(Call("Type: error, Time: " + new Date().toString() + ", NameMicroservice: " + name_microservice + ", NameClass: " + name_class + ", Event: Error in " + event + ", Error: " + TextError )) {
            System.out.println("Error");
            return gson.toJson(new SQLResult("Error"));
        } else {
            System.out.println("Error log");
            return gson.toJson(new SQLResult("ErrorLog"));
        }
    }

    public String getErrorHttp() {
        return ErrorHttp;
    }

    public void setErrorHttp(String ErrorHttp) {
        this.ErrorHttp = ErrorHttp;
    }
}
