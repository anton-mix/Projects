package Managers;

import ImageObjects.LoadImage;
import ImageObjects.LoadImageConverter;
import SQLClass.HttpSettings.HttpSettings;
import SQLClass.HttpSettings.UserEmailsConverter;
import SQLClass.HttpSettings.UserEmails;
import SQLClass.HttpSettings.HttpSettingsConverter;
import SQLClass.GlobalVariables;
import MyEnum.RestMailingEnum;
import MyInterface.InterfaceConverter;

import RestObject.RestObject;
import RestObject.RestObjectConverter;
import SQLClass.Contact.*;
import SQLClass.Message.*;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.User.UserSystemInfo;
import SQLClass.User.UserSystemInfoConverter;
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
    String name_class = "MailingAnswerManager";
    
    ContactsManager contact_manager;
    GroupsContactsManager groupsContacts_manager;
    TypesGroupsContactsManager typesGroupsContacts_manager;
    MessagesManager messages_manager;
    HttpSettingsManager settings_manager;
    RamManager ram_manager;
    ImageManager image_manager;
    
    public AnswerManager(Channel channel, GlobalVariables global_variables, LogManager log_manager) {
        super(channel);
        this.channel = channel;
        this.global_variables = global_variables;
        this.log_manager = log_manager;
   
        try {
           
            Class.forName("com.mysql.jdbc.Driver").newInstance();
            connection_string = "jdbc:mysql://" + global_variables.getDb_ip_sql()+ "/" + global_variables.getDb_name()+ "?useSSL=false&user=" + global_variables.getDb_user_name() + "&password=" + global_variables.getDb_password() + "&useUnicode=true&characterEncoding=UTF-8&serverTimezone=UTC";             
            
            contact_manager = new ContactsManager(global_variables, null, log_manager);
            groupsContacts_manager = new GroupsContactsManager(global_variables, null, log_manager);
            typesGroupsContacts_manager = new TypesGroupsContactsManager(global_variables, null, log_manager);
            messages_manager = new MessagesManager(global_variables, null, log_manager);
            settings_manager = new HttpSettingsManager(global_variables, null, log_manager);
            message_manager = new MessageManager(global_variables);
            image_manager = new ImageManager(log_manager, global_variables);
            
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(ContactModel.class, new ContactModelConverter());
            builder.registerTypeAdapter(ViewContact.class, new ViewContactConverter());
            builder.registerTypeAdapter(FilterContact.class, new FilterContactConverter());
            builder.registerTypeAdapter(ViewGroupContacts.class, new ViewGroupContactsConverter());
            builder.registerTypeAdapter(GroupContacts.class, new GroupContactsConverter());
            builder.registerTypeAdapter(TypeGroupContacts.class, new TypeGroupContactsConverter());
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());  
            builder.registerTypeAdapter(UserSystemInfo.class, new UserSystemInfoConverter());
            builder.registerTypeAdapter(FilterMessage.class, new FilterMessageConverter());
            builder.registerTypeAdapter(HttpSettings.class, new HttpSettingsConverter());
            builder.registerTypeAdapter(UserEmails.class, new UserEmailsConverter());
            builder.registerTypeAdapter(LoadImage.class, new LoadImageConverter());
            builder.registerTypeAdapter(SingleMessagesModel.class, new SingleMessagesModelConverter());
            builder.registerTypeAdapter(ViewSingleMessage.class, new ViewSingleMessageConverter());
            gson = builder.create();       
           
            ram_manager = new RamManager(log_manager, message_manager, gson);
        } 
        catch (Exception ex) {
            log_manager.CallError(ex.toString(), "AnswerManagerMailing. Error on init", name_class);
        } 
    }

    @Override
    public void handleDelivery(String string, Envelope envlp, AMQP.BasicProperties bp, byte[] bytes) throws IOException {     
        AMQP.BasicProperties replyProps = new AMQP.BasicProperties.Builder()
                .correlationId(bp.getCorrelationId())
                .build();
        ResData response = new ResData();
        response.setRes(new String(bytes, "UTF-8"));
        
        SQLResult res;

        try {
            RestObject rest_obj = gson.fromJson(response.getRes(), RestObject.class);
            
            log_manager.CallEvent(" receive message: " + rest_obj.toString(RestMailingEnum.values()[rest_obj.getCode()].name()), name_class);
            System.out.println("receive message: " + rest_obj.toString(RestMailingEnum.values()[rest_obj.getCode()].name()));
            
            connection_with_db = DriverManager.getConnection(connection_string);
            connection_with_db.setAutoCommit(false);
            contact_manager.setConnection_with_db(connection_with_db);
            groupsContacts_manager.setConnection_with_db(connection_with_db);
            typesGroupsContacts_manager.setConnection_with_db(connection_with_db);
            messages_manager.setConnection_with_db(connection_with_db);
            settings_manager.setConnection_with_db(connection_with_db);
            
            contact_manager.setErrorHttp(RestMailingEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            groupsContacts_manager.setErrorHttp(RestMailingEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            typesGroupsContacts_manager.setErrorHttp(RestMailingEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
            messages_manager.setErrorHttp(RestMailingEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            log_manager.setErrorHttp(RestMailingEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
            settings_manager.setErrorHttp(RestMailingEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
            image_manager.setErrorHttp(RestMailingEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
       
            res = ram_manager.CheckSessionKey(rest_obj.getSession_key());
            if (res.getId() == -1 || !"Ok".equals(res.getResult())) {
                switch (RestMailingEnum.values()[rest_obj.getCode()]) {
                    case Http_settings_get_path_image:
                        response = image_manager.getPathImage(gson.fromJson(rest_obj.getJson_message(), LoadImage.class));
                        break;
                    default:
                        response = new ResData(false, "{\"result\":\"auth-failed\"}");
                        break;
                }
            } else {
                boolean need_commit = false;
                rest_obj.getJson_message().addProperty("user_id", res.getId());
                switch (RestMailingEnum.values()[rest_obj.getCode()]) {
                    case Contact_select_one:
                        response = contact_manager.getContact(gson.fromJson(rest_obj.getJson_message(), ViewContact.class));
                        break; 
                    case Contact_select_view:
                        response = contact_manager.getViewContacts(gson.fromJson(rest_obj.getJson_message(), FilterContact.class));
                        break;
                    case Contact_select_count:
                        response = contact_manager.getCountContacts(gson.fromJson(rest_obj.getJson_message(), FilterContact.class));
                        break; 
                    case Contact_select_with_count: 
                        response = contact_manager.Contact_select_with_count(gson.fromJson(rest_obj.getJson_message(), FilterContact.class));
                        break;
                    case Contact_add:
                        response = contact_manager.AddContact(gson.fromJson(rest_obj.getJson_message(), ContactModel.class));
                        need_commit = true;
                        break;     
                    case Contact_update:
                        response = contact_manager.UpdateContact(gson.fromJson(rest_obj.getJson_message(), ContactModel.class));
                        need_commit = true;
                        break; 
                    case Contact_update_status:
                        response = contact_manager.UpdateStatusContact(gson.fromJson(rest_obj.getJson_message(), ContactModel.class));
                        need_commit = true;
                        break;
                    case Contact_delete:
                        response = contact_manager.DeleteContact(gson.fromJson(rest_obj.getJson_message(), ViewContact.class));
                        need_commit = true;
                        break;
                    case Group_contact_select_view:
                        response = groupsContacts_manager.getGroupsContactsView(gson.fromJson(rest_obj.getJson_message(), ViewGroupContacts.class));
                        break;
                    case Group_contact_select_one:
                        response = groupsContacts_manager.getGroupContacts(gson.fromJson(rest_obj.getJson_message(), ViewGroupContacts.class));
                        break;
                    case Group_contact_add:
                        response = groupsContacts_manager.AddGroupsContacts(gson.fromJson(rest_obj.getJson_message(), GroupContacts.class));
                        need_commit = true;
                        break;
                    case Group_contact_delete: 
                        response = groupsContacts_manager.DeleteGroupsContacts(gson.fromJson(rest_obj.getJson_message(), GroupContacts.class));
                        need_commit = true;
                        break;
                    case Group_contact_update: 
                        response = groupsContacts_manager.UpdateGroupsContacts(gson.fromJson(rest_obj.getJson_message(), GroupContacts.class));
                        need_commit = true;
                        break;
                    case Types_groups_contacts_select:
                        response = typesGroupsContacts_manager.getTypesGroupsContacts(gson.fromJson(rest_obj.getJson_message(), TypeGroupContacts.class));
                        break;
                    case Types_groups_contacts_add:
                        response = typesGroupsContacts_manager.AddTypesGroupsContacts(gson.fromJson(rest_obj.getJson_message(), TypeGroupContacts.class));
                        need_commit = true;
                        break;
                    case Types_groups_contacts_delete: 
                        response = typesGroupsContacts_manager.DeleteTypesGroupsContacts(gson.fromJson(rest_obj.getJson_message(), TypeGroupContacts.class));   
                        need_commit = true;
                        break;
                    case Types_groups_contacts_update:
                        response = typesGroupsContacts_manager.UpdateTypesGroupsContacts(gson.fromJson(rest_obj.getJson_message(), TypeGroupContacts.class));
                        need_commit = true;
                        break;
                    case Messages_select_one:
                        response = messages_manager.getMessages(gson.fromJson(rest_obj.getJson_message(), SingleMessagesModel.class));
                        break;
                    case Messages_select_view:
                        response = messages_manager.getViewMessages(gson.fromJson(rest_obj.getJson_message(), FilterMessage.class));
                        break;
                    case Messages_select_count:
                        response = messages_manager.getCountMessages(gson.fromJson(rest_obj.getJson_message(), FilterMessage.class));
                        break;
                    case Messages_select_with_count: 
                        response = messages_manager.Messages_select_with_count(gson.fromJson(rest_obj.getJson_message(), FilterMessage.class));
                        break;
                    case Message_add:
                        response = messages_manager.AddMessage(gson.fromJson(rest_obj.getJson_message(), SingleMessagesModel.class));
                        need_commit = true;
                        break;
                    case Message_delete:
                        response = messages_manager.DeleteMessage(gson.fromJson(rest_obj.getJson_message(), ViewSingleMessage.class));
                        need_commit = true;
                        break;
                    case Message_update:
                        response = messages_manager.UpdateMessage(gson.fromJson(rest_obj.getJson_message(), SingleMessagesModel.class));
                        need_commit = true;
                        break;
                    case Http_settings_select:
                        response = settings_manager.getHttpSettings(gson.fromJson(rest_obj.getJson_message(), HttpSettings.class));
                        break;
                    case Http_settings_add:
                        response = settings_manager.AddHttpSettings(gson.fromJson(rest_obj.getJson_message(), HttpSettings.class));
                        need_commit = true;
                        break;
                    case Http_settings_delete:
                        response = settings_manager.DeleteHttpSettings(gson.fromJson(rest_obj.getJson_message(), HttpSettings.class));
                        need_commit = true;
                        break;
                    case Http_settings_update:
                        response = settings_manager.UpdateHttpSettings(gson.fromJson(rest_obj.getJson_message(), HttpSettings.class));
                        need_commit = true;
                        break;
                    case Http_settings_set_path_image:
                        response = image_manager.setPathImage(gson.fromJson(rest_obj.getJson_message(), LoadImage.class));
                        break;
                }
                if (need_commit) {
                    if (response.getIs_success()) {
                        connection_with_db.commit();
                    } else {
                        connection_with_db.rollback();
                    }
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
