package Managers;

import RestObject.RestObject;
import SQLClass.GlobalVariables;
import SQLClass.ResData;
import com.google.gson.Gson;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.rabbitmq.client.ConnectionFactory;
import com.rabbitmq.client.Connection;
import com.rabbitmq.client.Channel;
import com.rabbitmq.client.DefaultConsumer;
import com.rabbitmq.client.AMQP;
import com.rabbitmq.client.Envelope;

import java.io.IOException;
import java.util.UUID;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

public class MessageManager {
    private Connection connection;
    private Channel channel;
    GlobalVariables global_variables;
    
    public MessageManager(GlobalVariables global_variables) throws IOException, TimeoutException {
        this.global_variables = global_variables;
        ConnectionFactory factory = new ConnectionFactory();
        factory.setHost(global_variables.getHost());
        factory.setPort(global_variables.getPort());
        connection = factory.newConnection();
        channel = connection.createChannel(); 
    }
    
    public ResData call(Object message, Gson gson, String SessionKey, int Code,  String requestQueueName) throws IOException, InterruptedException {
        JsonParser parser = new JsonParser();
        JsonObject json_message = parser.parse(gson.toJson(message)).getAsJsonObject();
        return call(gson.toJson(new RestObject(json_message, SessionKey, Code), RestObject.class),requestQueueName);
    }
    
    public ResData call(String message, Gson gson, String SessionKey, int Code,  String requestQueueName) throws IOException, InterruptedException {
        JsonParser parser = new JsonParser();
        JsonObject json_message = parser.parse(message).getAsJsonObject();
        return call(gson.toJson(new RestObject(json_message, SessionKey, Code), RestObject.class),requestQueueName);
    }
    
    public ResData call(JsonObject message, Gson gson, String SessionKey, int Code, String CodeName,  String requestQueueName) throws IOException, InterruptedException {
        return call(gson.toJson(new RestObject(message, SessionKey, Code), RestObject.class),requestQueueName);
    }
    
    public ResData call(String message, String requestQueueName) throws IOException, InterruptedException {
        final String corrId = UUID.randomUUID().toString();
        String replyQueueName = channel.queueDeclare().getQueue();
        AMQP.BasicProperties props = new AMQP.BasicProperties
            .Builder()
            .correlationId(corrId)
            .replyTo(replyQueueName)
            .build();
        channel.basicPublish("", requestQueueName, props, message.getBytes("UTF-8"));

        final BlockingQueue<String> response = new ArrayBlockingQueue<String>(1);
        String ctag = channel.basicConsume(replyQueueName, true, new DefaultConsumer(channel) {
            @Override
            public void handleDelivery(String consumerTag, Envelope envelope, AMQP.BasicProperties properties, byte[] body) throws IOException {
                if (properties.getCorrelationId().equals(corrId)) {
                    response.offer(new String(body, "UTF-8"));
                }
            }
        });

        String result = response.take();
        /*
        String result = response.poll(10000, TimeUnit.MILLISECONDS);
        if(result == null)
            result = "{\"result\":\"Internal Error\"}";
        */
        channel.basicCancel(ctag);
        return new ResData(true, result);
    }

    public void close() throws IOException {
        connection.close();
    } 
}
