package Managers;

import SQLClass.GlobalVariables;
import SQLClass.GroupProducts;
import SQLClass.InvoicesToProductsItem;
import SQLClass.Product;
import SQLClass.ProductsToGroupProductsItem;
import com.rabbitmq.client.AMQP;
import com.rabbitmq.client.Channel;
import com.rabbitmq.client.DefaultConsumer;
import com.rabbitmq.client.Envelope;
import java.io.IOException;
import java.sql.SQLException;
import java.util.logging.Level;
import java.util.logging.Logger;


public class AnswerManager extends DefaultConsumer {
    LogManager log_manager;
    ProductsManager products_manager;
    GroupsProductsManager groups_products_manager;
    InvoicesToProductsManager invoices_to_products_manager;
    ProductsToGroupProductsManager products_to_group_products_manager;
    GlobalVariables global_variables;
    Channel channel;
    
    public AnswerManager(Channel channel, GlobalVariables global_variables ) {
        super(channel);
        
        this.channel = channel;
        this.global_variables = global_variables;
        
        log_manager = new LogManager(global_variables.getHost(), global_variables.getQueueNameLog(), "AnswerManagerProducts");
        
        try {
            products_manager = new ProductsManager(global_variables);
            groups_products_manager = new GroupsProductsManager(global_variables);
            invoices_to_products_manager = new InvoicesToProductsManager(global_variables);
            products_to_group_products_manager = new ProductsToGroupProductsManager(global_variables);
        } catch (SQLException ex) {
            Logger.getLogger(AnswerManager.class.getName()).log(Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            Logger.getLogger(AnswerManager.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            Logger.getLogger(AnswerManager.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ClassNotFoundException ex) {
            Logger.getLogger(AnswerManager.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    @Override
    public void handleDelivery(String string, Envelope envlp, AMQP.BasicProperties bp, byte[] bytes) throws IOException {        
        log_manager.CallEvent("receive from client");

        AMQP.BasicProperties replyProps = new AMQP.BasicProperties.Builder()
                .correlationId(bp.getCorrelationId())
                .build();

        String string_json = "";
        String response = "";
        try {
            response = new String(bytes, "UTF-8");
            String[] part_string = response.split(" ");
            
            System.out.println(response);
            

            
            switch (Integer.valueOf(part_string[0])) {
                case 0:
                    response = products_manager.getProducts();
                    break;
                case 1:
                    for (int i = 1; i < part_string.length - 2; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
                    Product pr = new Product(string_json);
                    response = products_manager.getProducts(new Product(string_json), Integer.valueOf(part_string[part_string.length - 2]), Integer.valueOf(part_string[part_string.length - 1]));
                    break;
                case 2:
                    for (int i = 1; i < part_string.length; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
 
                    response = products_manager.getCountProducts(new Product(string_json));

                    break;
                case 3:
                    for (int i = 1; i < part_string.length; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
                    
                    response = products_manager.AddProduct(new Product(string_json));
                    
                    break;
                case 4:
                    
                    response = products_manager.DeleteProduct(Integer.valueOf(part_string[1]));
                    
                    break;
                case 5:
                    for (int i = 1; i < part_string.length; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
                    
                    response = products_manager.UpdateProduct(new Product(string_json));
                    
                    break;
                case 6:
                    for (int i = 1; i < part_string.length - 2; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
                    
                    response = groups_products_manager.getGroupsProducts(new GroupProducts(string_json), Integer.valueOf(part_string[part_string.length - 2]), Integer.valueOf(part_string[part_string.length - 1]));
                    
                    break;
                case 7:
                    for (int i = 1; i < part_string.length; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
 
                    response = groups_products_manager.getCountGroupsProducts(new GroupProducts(string_json));

                    break;
                case 8:
                    for (int i = 1; i < part_string.length; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
                    
                    response = groups_products_manager.AddGroupProducts(new GroupProducts(string_json));
                    
                    break;
                case 9:
                    response = groups_products_manager.DeleteGroupProducts(Integer.valueOf(part_string[1]));
                    break;
                case 10:
                    for (int i = 1; i < part_string.length; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
                    
                    response = groups_products_manager.UpdateGroupProducts(new GroupProducts(string_json));
                break;
                case 11:
                    for (int i = 1; i < part_string.length - 2; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
                    
                    response = invoices_to_products_manager.getInvoicesToProductsItems(new InvoicesToProductsItem(string_json), Integer.valueOf(part_string[part_string.length - 2]), Integer.valueOf(part_string[part_string.length - 1]));
                    
                    break;
                case 12:
                    for (int i = 1; i < part_string.length; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
 
                    response = invoices_to_products_manager.getCountInvoicesToProductsItems(new InvoicesToProductsItem(string_json));

                    break;
                case 13:
                    for (int i = 1; i < part_string.length; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
                    
                    response = invoices_to_products_manager.AddInvoicesToProductsItem(new InvoicesToProductsItem(string_json));
                    
                    break;
                case 14:
                    
                    response = invoices_to_products_manager.DeleteInvoicesToProductsItem(Integer.valueOf(part_string[1]));
                    
                    break;
                case 15:
                    for (int i = 1; i < part_string.length; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
                    
                    response = invoices_to_products_manager.UpdateInvoicesToProductsItem(new InvoicesToProductsItem(string_json));
                break;
                case 16:
                    for (int i = 1; i < part_string.length - 2; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
                    
                    response = products_to_group_products_manager.getProductsToGroupProductsItems(new ProductsToGroupProductsItem(string_json), Integer.valueOf(part_string[part_string.length - 2]), Integer.valueOf(part_string[part_string.length - 1]));
                    
                    break;
                case 17:
                    for (int i = 1; i < part_string.length; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
 
                    response = products_to_group_products_manager.getCountProductsToGroupProductsItems(new ProductsToGroupProductsItem(string_json));

                    break;
                case 18:
                    for (int i = 1; i < part_string.length; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
                    
                    response = products_to_group_products_manager.AddProductsToGroupProductsItem(new ProductsToGroupProductsItem(string_json));
                    
                    break;
                case 19:
 
                    response = products_to_group_products_manager.DeleteProductsToGroupProductsItem(Integer.valueOf(part_string[1]));
                    
                    break;
                case 20:
                    for (int i = 1; i < part_string.length; i++) {
                        string_json = string_json + " " + part_string[i];
                    }
                    
                    response = products_to_group_products_manager.UpdateProductsToGroupProductsItem(new ProductsToGroupProductsItem(string_json));
                break;
                case 21:
                    response = products_manager.getProductsByGroupId(Integer.valueOf(part_string[1]));
                break;
            }
            System.out.println(response);
            
        } catch (Exception ex) {
            log_manager.CallError(ex.toString(), "call to queue");
        } finally {
            log_manager.CallEvent("send to client");

            channel.basicPublish("", bp.getReplyTo(), replyProps, response.getBytes("UTF-8"));
            channel.basicAck(envlp.getDeliveryTag(), false);
            synchronized (this) {
                this.notify();
            }
        }
    }
}
