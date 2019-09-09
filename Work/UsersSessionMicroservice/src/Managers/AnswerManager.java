package Managers;

import ImageObjects.LoadImage;
import ImageObjects.LoadImageConverter;
import MyInterface.InterfaceConverter;
import MessageClass.*;
import MyEnum.*;
import SQLClass.User.*;
import SQLClass.*;
import RestObject.*;
import QueueClass.*;
import MyException.MyNotAuthorizedException;
import SQLClass.Contact.TypeGroupContacts;
import SQLClass.Contact.TypeGroupContactsConverter;
import SQLClass.HttpSettings.UserEmails;
import SQLClass.HttpSettings.UserEmailsConverter;
import SQLClass.Product.GroupProducts;
import SQLClass.Product.GroupProductsConverter;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.rabbitmq.client.AMQP;
import com.rabbitmq.client.Channel;
import com.rabbitmq.client.DefaultConsumer;
import com.rabbitmq.client.Envelope;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

public class AnswerManager extends DefaultConsumer implements InterfaceConverter {
    private Channel channel;
    Connection connection_with_db = null;
    GlobalVariables global_variables;
    LogManager log_manager;
    MessageManager message_manager;
    Gson gson;
    
    String connection_string;
    String name_class = "UserAnswerManager";
    boolean need_commit = false;
    
    UsersManager users_manager;
    UserSystemInfoManager userSystemInfo_manager;
    RamManager ram_manager;
    ImageManager loadImage_manager;

    public AnswerManager(Channel channel, GlobalVariables global_variables, LogManager log_manager) {
        super(channel);
        this.channel = channel;
        this.global_variables = global_variables;
        this.log_manager = log_manager;
        
        try {
            Class.forName("com.mysql.jdbc.Driver").newInstance();
            connection_string = "jdbc:mysql://" + global_variables.getDb_ip_sql() + "/" + global_variables.getDb_name() + "?useSSL=false&user=" + global_variables.getDb_user_name() + "&password=" + global_variables.getDb_password() + "&useUnicode=true&characterEncoding=UTF-8&serverTimezone=UTC";             
            System.out.println(connection_string);
            users_manager = new UsersManager(global_variables, null, log_manager);
            userSystemInfo_manager = new UserSystemInfoManager(global_variables, null, log_manager);
            message_manager = new MessageManager(global_variables);
            loadImage_manager = new ImageManager(log_manager, global_variables);
            
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(UserModel.class, new UserModelConverter());
            builder.registerTypeAdapter(UserSystemInfo.class, new UserSystemInfoConverter());
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());  
            builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());  
            builder.registerTypeAdapter(RessiveClass.class, new RessiveClassConverter()); 
            builder.registerTypeAdapter(MessageAuth.class, new MessageAuthConverter()); 
            builder.registerTypeAdapter(UserModel.class, new UserModelConverter()); 
            builder.registerTypeAdapter(GroupProducts.class, new GroupProductsConverter());
            builder.registerTypeAdapter(TypeGroupContacts.class, new TypeGroupContactsConverter());
            builder.registerTypeAdapter(MessageChangePassword.class, new MessageChangePasswordConverter());
            builder.registerTypeAdapter(UserUpdatePassword.class, new UserUpdatePasswordConverter());
            builder.registerTypeAdapter(ArrUserPasswordKeyValue.class, new ArrUserPasswordKeyValueConverter());
            builder.registerTypeAdapter(UserEmails.class, new UserEmailsConverter());
            builder.registerTypeAdapter(LoadImage.class, new LoadImageConverter());
            gson = builder.create();    
            
            ram_manager = new RamManager(log_manager, message_manager, gson);
        } catch (Exception ex) {
            log_manager.CallError(ex.toString(), "AnswerManager. Error on init", name_class);
        } 
    }
    
    @Override
    public void handleDelivery(String string, Envelope envlp, AMQP.BasicProperties bp, byte[] bytes) throws IOException {     
        AMQP.BasicProperties replyProps = new AMQP.BasicProperties.Builder()
            .correlationId(bp.getCorrelationId())
            .build();
        
        ResData response = new ResData();
        response.setRes(new String(bytes, "UTF-8"));
        SQLResult resQueue;
        
        log_manager.setErrorHttp(response.getRes());
       
        
        try {
            RestObject rest_obj = gson.fromJson(response.getRes(), RestObject.class);
            System.out.println("receive message: " + rest_obj.toString(RestUserEnum.values()[rest_obj.getCode()].name()));
            log_manager.CallEvent(" receive message: " + rest_obj.toString(RestUserEnum.values()[rest_obj.getCode()].name()), name_class);
            
            connection_with_db = DriverManager.getConnection(connection_string);
            connection_with_db.setAutoCommit(false);
            users_manager.setConnection_with_db(connection_with_db);
            userSystemInfo_manager.setConnection_with_db(connection_with_db);
            
            users_manager.setErrorHttp(RestUserEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            userSystemInfo_manager.setErrorHttp(RestUserEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            log_manager.setErrorHttp(RestUserEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            loadImage_manager.setErrorHttp(RestUserEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());

            resQueue = ram_manager.CheckSessionKey(rest_obj.getSession_key());
            
            if (resQueue.getId() == -1 || !"Ok".equals(resQueue.getResult())) {
                switch (RestUserEnum.values()[rest_obj.getCode()]) {
                    case User_select_one:
                        response = users_manager.getOneUser(gson.fromJson(rest_obj.getJson_message(), UserModel.class), ram_manager);
                        break;
                    case User_add:
                        response = users_manager.AddUser(gson.fromJson(rest_obj.getJson_message(), UserAuthModel.class), ram_manager, userSystemInfo_manager);
                        need_commit = true;
                        break;
                    case User_pre_reg:
                        response = users_manager.UserPreReg(gson.fromJson(rest_obj.getJson_message(), UserAuthModel.class), ram_manager, message_manager);
                        break;
                    case User_check_param:
                        response = users_manager.UserCheckParam(gson.fromJson(rest_obj.getJson_message(), UserModel.class));
                        break;
                    case User_begin_change_password:
                        response = users_manager.UserBeginChangePassword(gson.fromJson(rest_obj.getJson_message(), UserModel.class), ram_manager, message_manager);
                        break;
                    case User_end_change_password:
                        response = users_manager.UserEndChangePassword(gson.fromJson(rest_obj.getJson_message(), UserModel.class), ram_manager);
                        need_commit = true;
                        break;
                    case User_get_path_image:
                        response = loadImage_manager.getPathImage(gson.fromJson(rest_obj.getJson_message(), LoadImage.class));
                        break;  
                    default:
                        response = new ResData(false, "{\"result\":\"auth-failed\"}");
                        break;
                }
            } else {
                rest_obj.getJson_message().addProperty("user_id", resQueue.getUser_id());
                rest_obj.getJson_message().addProperty("id", resQueue.getId());
                switch (RestUserEnum.values()[rest_obj.getCode()]) {
                    case User_select_all:
                        response = users_manager.getUser(gson.fromJson(rest_obj.getJson_message(), UserModel.class));
                        break;
                    case User_cancel:
                        response = users_manager.UserCancel(gson.fromJson(rest_obj.getJson_message(), UserAuthModel.class), ram_manager);
                        break;
                    case User_delete:
                        response = users_manager.DeleteUser(gson.fromJson(rest_obj.getJson_message(), UserModel.class));
                        need_commit = true;
                        break;
                    case User_update:
                        response = users_manager.UpdateUser(gson.fromJson(rest_obj.getJson_message(), UserModel.class));
                        need_commit = true;
                        break;
                    case User_select_emails:
                        response = users_manager.getUserEmails(gson.fromJson(rest_obj.getJson_message(), UserEmails.class));
                        break;
                    case User_update_account_password:
                        response = users_manager.UserUpdateAccountPassword(gson.fromJson(rest_obj.getJson_message(), UserUpdatePassword.class));
                        need_commit = true;
                        break;
                    case User_set_path_image:
                        response = loadImage_manager.setPathImage(gson.fromJson(rest_obj.getJson_message(), LoadImage.class));
                        break;
                }
            }
            if (need_commit) {
                if (response.getIs_success()) {
                    connection_with_db.commit();
                } else {
                    connection_with_db.rollback();
                }
            }
        }
        catch(Exception e) {
            response = log_manager.CallError(e.toString(), "handleDelivery(); Error to send.", name_class);
        }
        finally {
            try {
                connection_with_db.close();
            } catch (SQLException ex) {
                log_manager.CallEvent("handleDelivery(); Connection Close Error", name_class);
            }
            catch (Exception ex) {
                log_manager.CallEvent("handleDelivery(); Connection Close Error", name_class);
            }
            if(bp.getReplyTo() == null) {
                log_manager.CallEvent("send message: " + response.getRes() + " is drop", name_class);
                System.out.println("send message: " + response.getRes() + " is drop");
                channel.basicAck(envlp.getDeliveryTag(), false);
            }
            else {
                log_manager.CallEvent("send message: " + response.getRes(), name_class);
                System.out.println("send message: " + response.getRes());
                channel.basicPublish("", bp.getReplyTo(), replyProps, response.getRes().getBytes("UTF-8"));
                channel.basicAck(envlp.getDeliveryTag(), false);
            }
            synchronized (this) {
                this.notify();
            }
        }
    } 
}
