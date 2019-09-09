package Managers;

import ImageObjects.LoadImage;
import ImageObjects.LoadImageConverter;
import MyEnum.RestShopEnum;
import SQLClass.GlobalVariables;
import SQLClass.Invoice.*;
import MyInterface.InterfaceConverter;

import RestObject.RestObject;
import RestObject.RestObjectConverter;
import SQLClass.Product.*;
import SQLClass.ReminderMaill.ReminderMessage;
import SQLClass.ReminderMaill.ReminderMessageConverter;
import SQLClass.ReminderPayments.ReminderPaymentsModel;
import SQLClass.ReminderPayments.ReminderPaymentsModelConverter;
import SQLClass.ReminderSurcharge.ReminderSurchargeModel;
import SQLClass.ReminderSurcharge.ReminderSurchargeModelConverter;
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
    String name_class = "ShopAnswerManager";
    
    InvoicesManager invoices_manager;
    InvoiceToProductsManager invoicesToProducts_manager;
    InvoiceEventManager invoiceEvent_manager;
    ProductsManager products_manager;
    GroupsProductsManager groupsProducts_manager;
    RamManager ram_manager;
    ImageManager image_manager;
    ReminderPaymentsManager reminderPayments_manager;
    ReminderPaymentsMessageManager reminderPaymentsMessage_manager;
    ReminderSurchargeManager reminderSurcharge_manager;
    ReminderSurchargeMessageManager reminderSurchargeMessage_manager;
    
    public AnswerManager(Channel channel, GlobalVariables global_variables, LogManager log_manager) {
        super(channel);
        this.channel = channel;
        this.global_variables = global_variables;
        this.log_manager = log_manager;
   
        try {
            Class.forName("com.mysql.jdbc.Driver").newInstance();
            connection_string = "jdbc:mysql://" + global_variables.getDb_ip_sql() + "/" + global_variables.getDb_name() + "?useSSL=false&user=" + global_variables.getDb_user_name() + "&password=" + global_variables.getDb_password() + "&useUnicode=true&characterEncoding=UTF-8&serverTimezone=UTC";             
            
            message_manager = new MessageManager(global_variables);
            invoices_manager = new InvoicesManager(global_variables, null, log_manager);
            invoicesToProducts_manager = new InvoiceToProductsManager(global_variables, null, log_manager);
            invoiceEvent_manager = new InvoiceEventManager(global_variables, null, log_manager);
            products_manager = new ProductsManager(global_variables, null, log_manager);
            groupsProducts_manager = new GroupsProductsManager(global_variables, null, log_manager);
            image_manager = new ImageManager(log_manager, global_variables, products_manager);
            reminderPayments_manager = new ReminderPaymentsManager(global_variables, connection_with_db, log_manager);
            reminderSurcharge_manager = new ReminderSurchargeManager(global_variables, connection_with_db, log_manager);
            reminderPaymentsMessage_manager = new ReminderPaymentsMessageManager(global_variables, connection_with_db, log_manager);
            reminderSurchargeMessage_manager = new ReminderSurchargeMessageManager(global_variables, connection_with_db, log_manager);
            
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(InvoiceModel.class, new InvoiceModelConverter());
            builder.registerTypeAdapter(ViewInvoice.class, new ViewInvoiceConverter());
            builder.registerTypeAdapter(InvoiceEvent.class, new InvoiceEventConverter());
            builder.registerTypeAdapter(ProductModel.class, new ProductModelConverter());
            builder.registerTypeAdapter(ViewProduct.class, new ViewProductConverter());
            builder.registerTypeAdapter(FilterProduct.class, new FilterProductConverter());
            builder.registerTypeAdapter(GroupProducts.class, new GroupProductsConverter());
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter()); 
            builder.registerTypeAdapter(UserSystemInfo.class, new UserSystemInfoConverter());
            builder.registerTypeAdapter(LoadImage.class, new LoadImageConverter());
            builder.registerTypeAdapter(ReminderPaymentsModel.class, new ReminderPaymentsModelConverter());
            builder.registerTypeAdapter(ReminderSurchargeModel.class, new ReminderSurchargeModelConverter());
            builder.registerTypeAdapter(ReminderMessage.class, new ReminderMessageConverter());
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
        SQLResult res;
            
        log_manager.setErrorHttp(response.getRes());
        
        try {
            RestObject rest_obj = gson.fromJson(response.getRes(), RestObject.class);
            
            System.out.println("receive message: " + rest_obj.toString(RestShopEnum.values()[rest_obj.getCode()].name()));
            log_manager.CallEvent(" receive message: " + rest_obj.toString(RestShopEnum.values()[rest_obj.getCode()].name()), name_class);
            
            connection_with_db = DriverManager.getConnection(connection_string);
            connection_with_db.setAutoCommit(false);
            
            invoices_manager.setConnection_with_db(connection_with_db);
            invoicesToProducts_manager.setConnection_with_db(connection_with_db);
            invoiceEvent_manager.setConnection_with_db(connection_with_db);
            products_manager.setConnection_with_db(connection_with_db);
            groupsProducts_manager.setConnection_with_db(connection_with_db);
            reminderPayments_manager.setConnection_with_db(connection_with_db);
            reminderSurcharge_manager.setConnection_with_db(connection_with_db);
            reminderPaymentsMessage_manager.setConnection_with_db(connection_with_db);
            reminderSurchargeMessage_manager.setConnection_with_db(connection_with_db);
            
            
            invoices_manager.setErrorHttp(RestShopEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            invoicesToProducts_manager.setErrorHttp(RestShopEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            invoiceEvent_manager.setErrorHttp(RestShopEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            products_manager.setErrorHttp(RestShopEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            groupsProducts_manager.setErrorHttp(RestShopEnum.values()[rest_obj.getCode()].name()  + " " + rest_obj.getJson_message().toString());
            log_manager.setErrorHttp(RestShopEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
            image_manager.setErrorHttp(RestShopEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
            reminderPayments_manager.setErrorHttp(RestShopEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
            reminderSurcharge_manager.setErrorHttp(RestShopEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
            reminderPaymentsMessage_manager.setErrorHttp(RestShopEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
            reminderSurchargeMessage_manager.setErrorHttp(RestShopEnum.values()[rest_obj.getCode()].name() + " " + rest_obj.getJson_message().toString());
            
            res = ram_manager.CheckSessionKey(rest_obj.getSession_key());
            if (res.getId() == -1 || !"Ok".equals(res.getResult())) {
                switch (RestShopEnum.values()[rest_obj.getCode()]) {
                    case Product_get_path_image:
                        response = image_manager.getPathImage(gson.fromJson(rest_obj.getJson_message(), LoadImage.class));
                        break;
                    default:
                        response = new ResData(false, "{\"result\":\"auth-failed\"}");
                        break;
                }
            } 
            else {
                rest_obj.getJson_message().addProperty("user_id", res.getId());
                boolean need_commit = false;
                switch (RestShopEnum.values()[rest_obj.getCode()]) {
                    case Invoice_select_all: 
                        response = invoices_manager.getInvoices();
                        break;
                    case Invoice_select_one: 
                        response = invoices_manager.getInvoice(gson.fromJson(rest_obj.getJson_message(), ViewInvoice.class));
                        break;
                    case Invoice_select_view: 
                        response = invoices_manager.getViewInvoices(gson.fromJson(rest_obj.getJson_message(), FilterInvoice.class));
                        break;
                    case Invoice_select_count: 
                        response = invoices_manager.getCountInvoices(gson.fromJson(rest_obj.getJson_message(), FilterInvoice.class));
                        break;
                    case Invoice_select_with_count: 
                        response = invoices_manager.Invoice_select_with_count(gson.fromJson(rest_obj.getJson_message(), FilterInvoice.class));
                        break;
                    case Invoice_add:  
                        response = invoices_manager.AddInvoice(gson.fromJson(rest_obj.getJson_message(), InvoiceModel.class));
                        need_commit = true;
                        break;
                    case Invoice_delete:
                        response = invoices_manager.DeleteInvoice(gson.fromJson(rest_obj.getJson_message(), ViewInvoice.class));
                        need_commit = true;
                        break;
                    case Invoice_update:
                        response = invoices_manager.UpdateInvoice(gson.fromJson(rest_obj.getJson_message(), InvoiceModel.class));
                        need_commit = true;
                        break;
                        
                    case Invoice_event_select:
                        response = invoiceEvent_manager.getInvoiceEvents(gson.fromJson(rest_obj.getJson_message(), InvoiceEvent.class));
                        break;
                    case Invoice_event_select_count:
                        response = invoiceEvent_manager.getCountInvoiceEvents(gson.fromJson(rest_obj.getJson_message(), InvoiceEvent.class));
                        break;
                    case Invoice_event_add:
                        response = invoiceEvent_manager.AddInvoiceEvents(gson.fromJson(rest_obj.getJson_message(), InvoiceEvent.class));
                        need_commit = true;
                        break;
                    case Invoice_event_delete:
                        response = invoiceEvent_manager.DeleteInvoiceEvents(gson.fromJson(rest_obj.getJson_message(), InvoiceEvent.class));
                        need_commit = true;
                        break;
                    case Invoice_event_update:
                        response = invoiceEvent_manager.UpdateInvoiceEvents(gson.fromJson(rest_obj.getJson_message(), InvoiceEvent.class));
                        need_commit = true;
                        break;
                        
                    case Product_select_all:
                        response = products_manager.getProducts();
                        break;
                    case Product_select_count:
                        response = products_manager.getCountProducts(gson.fromJson(rest_obj.getJson_message(), FilterProduct.class));
                        break;
                    case Product_select_one:
                        response = products_manager.getProduct(gson.fromJson(rest_obj.getJson_message(), ViewProduct.class));
                        break;
                    case Product_select_view:
                        response = products_manager.getViewProducts(gson.fromJson(rest_obj.getJson_message(), FilterProduct.class));
                        break;
                    case Product_select_short:
                        response = products_manager.getShortProducts(gson.fromJson(rest_obj.getJson_message(), FilterProduct.class));
                        break;
                    case Product_select_with_count: 
                        response = products_manager.Product_select_with_count(gson.fromJson(rest_obj.getJson_message(), FilterProduct.class));
                        break;
                    case Product_add:
                        response = products_manager.AddProduct(gson.fromJson(rest_obj.getJson_message(), ProductModel.class));
                        need_commit = true;
                        break;
                    case Product_delete:
                        response = products_manager.DeleteProduct(gson.fromJson(rest_obj.getJson_message(), ProductModel.class));
                        need_commit = true;
                        break;
                    case Product_update:
                        response = products_manager.UpdateProduct(gson.fromJson(rest_obj.getJson_message(), ProductModel.class));
                        need_commit = true;
                        break;
                    case Product_set_path_image:
                        response = image_manager.setPathImage(gson.fromJson(rest_obj.getJson_message(), LoadImage.class));
                        break;
                        
                    case Group_product_select:
                        response = groupsProducts_manager.getGroupsProducts(gson.fromJson(rest_obj.getJson_message(), GroupProducts.class));
                        break;
                    case Group_product_select_count:
                        response = groupsProducts_manager.getCountGroupsProducts(gson.fromJson(rest_obj.getJson_message(), GroupProducts.class));
                        break;
                    case Group_product_add:
                        response = groupsProducts_manager.AddGroupsProducts(gson.fromJson(rest_obj.getJson_message(), GroupProducts.class));
                        need_commit = true;
                        break;
                    case Group_product_delete: 
                        response = groupsProducts_manager.DeleteGroupsProducts(gson.fromJson(rest_obj.getJson_message(), GroupProducts.class));
                        need_commit = true;
                        break;
                    case Group_product_update: 
                        response = groupsProducts_manager.UpdateGroupsProducts(gson.fromJson(rest_obj.getJson_message(), GroupProducts.class));
                        need_commit = true;
                        break;
                    case Reminder_payments_select:
                        response = reminderPayments_manager.getReminderPayments(gson.fromJson(rest_obj.getJson_message(), ReminderPaymentsModel.class));
                        break;
                    case Reminder_payments_select_count:
                        response = reminderPayments_manager.getCountReminderPayments(gson.fromJson(rest_obj.getJson_message(), ReminderPaymentsModel.class));
                        break;
                    case Reminder_payments_select_with_count:
                        response = reminderPayments_manager.getSelectWithCountReminderPayments(gson.fromJson(rest_obj.getJson_message(), ReminderPaymentsModel.class));
                        break;
                    case Reminder_payments_add:
                        response = reminderPayments_manager.AddReminderPayments(gson.fromJson(rest_obj.getJson_message(), ReminderPaymentsModel.class));
                        need_commit = true;
                        break;
                    case Reminder_payments_delete:
                        response = reminderPayments_manager.DeleteReminderPayments(gson.fromJson(rest_obj.getJson_message(), ReminderPaymentsModel.class));
                        need_commit = true;
                        break;
                    case Reminder_payments_update:
                        response = reminderPayments_manager.UpdateReminderPayments(gson.fromJson(rest_obj.getJson_message(), ReminderPaymentsModel.class));
                        need_commit = true;
                        break;
    
                    case Reminder_surcharge_select:
                        response = reminderSurcharge_manager.getCountReminderSurcharge(gson.fromJson(rest_obj.getJson_message(), ReminderSurchargeModel.class));
                        break;
                    case Reminder_surcharge_select_count:
                        response = reminderSurcharge_manager.getCountReminderSurcharge(gson.fromJson(rest_obj.getJson_message(), ReminderSurchargeModel.class));
                        break;
                    case Reminder_surcharge_select_with_count:
                        response = reminderSurcharge_manager.getSelectWithCountReminderPayments(gson.fromJson(rest_obj.getJson_message(), ReminderSurchargeModel.class));
                        break;
                    case Reminder_surcharge_add:
                        response = reminderSurcharge_manager.AddReminderSurcharge(gson.fromJson(rest_obj.getJson_message(), ReminderSurchargeModel.class));
                        need_commit = true;
                        break;
                    case Reminder_surcharge_delete:
                        response = reminderSurcharge_manager.DeleteReminderSurcharge(gson.fromJson(rest_obj.getJson_message(), ReminderSurchargeModel.class));
                        need_commit = true;
                        break;
                    case Reminder_surcharge_update:
                        response = reminderSurcharge_manager.UpdateReminderSurcharge(gson.fromJson(rest_obj.getJson_message(), ReminderSurchargeModel.class));
                        need_commit = true;
                        break;
                        
                    case Reminder_payments_messages_select:
                        reminderPaymentsMessage_manager.getReminderPaymentsMessage(gson.fromJson(rest_obj.getJson_message(), ReminderMessage.class));
                        break;
                    case Reminder_payments_messages_add:
                        reminderPaymentsMessage_manager.AddReminderPaymentsMessage(gson.fromJson(rest_obj.getJson_message(), ReminderMessage.class));
                        need_commit = true;
                        break;
                    case Reminder_payments_messages_delete:
                        reminderPaymentsMessage_manager.DeleteReminderPaymentsMessage(gson.fromJson(rest_obj.getJson_message(), ReminderMessage.class));
                        need_commit = true;
                        break;
                    case Reminder_payments_messages_update:
                        reminderPaymentsMessage_manager.UpdateReminderPaymentsMessage(gson.fromJson(rest_obj.getJson_message(), ReminderMessage.class));
                        need_commit = true;
                        break;
    
                    case Reminder_surcharge_messages_select:
                        reminderSurchargeMessage_manager.getReminderSurchargeMessage(gson.fromJson(rest_obj.getJson_message(), ReminderMessage.class));
                        break;
                    case Reminder_surcharge_messages_add:
                        reminderSurchargeMessage_manager.AddReminderSurchargeMessage(gson.fromJson(rest_obj.getJson_message(), ReminderMessage.class));
                        need_commit = true;
                        break;
                    case Reminder_surcharge_messages_delete:
                        reminderSurchargeMessage_manager.DeleteReminderSurchargeMessage(gson.fromJson(rest_obj.getJson_message(), ReminderMessage.class));
                        need_commit = true;
                        break;
                    case Reminder_surcharge_messages_update:
                        reminderSurchargeMessage_manager.UpdateReminderSurchargeMessage(gson.fromJson(rest_obj.getJson_message(), ReminderMessage.class));
                        need_commit = true;
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
            MyFinally(bp, response,  envlp, replyProps);
        }
        catch(Exception e) {
            response = log_manager.CallError(e.toString(), "handleDelivery(); Error to send.", name_class);
            MyFinally(bp, response,  envlp, replyProps);
        }
    }
    
    void MyFinally(AMQP.BasicProperties bp, ResData response, Envelope envlp, AMQP.BasicProperties replyProps) {
        try {
            connection_with_db.close();
        } catch (Exception ex) {
            log_manager.CallEvent("handleDelivery(); Connection Close Error", name_class);
        }

        try {
            if (bp.getReplyTo() == null) {
                log_manager.CallEvent("send message: " + response.getRes() + " is drop", name_class);
                System.out.println("send message: " + response.getRes() + " is drop");
                channel.basicAck(envlp.getDeliveryTag(), false);
            } else {
                log_manager.CallEvent("send message: " + response.getRes(), name_class);
                System.out.println("send message: " + response.getRes());
                channel.basicPublish("", bp.getReplyTo(), replyProps, response.getRes().getBytes("UTF-8"));
                channel.basicAck(envlp.getDeliveryTag(), false);
            }
        } catch (Exception ex) {
            log_manager.CallEvent("handleDelivery(); Connection Close Error", name_class);
        }
        
        synchronized (this) {
            this.notify();
        }
    }
}
