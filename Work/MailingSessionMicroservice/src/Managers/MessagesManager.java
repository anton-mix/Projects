
package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.GlobalVariables;
import SQLClass.Message.FilterMessage;
import SQLClass.Message.FilterMessageConverter;
import SQLClass.Message.MessageImageLinks;
import SQLClass.Message.MessageLinks;
import SQLClass.Message.SingleMessagesModel;
import SQLClass.Message.SingleMessagesModelConverter;
import SQLClass.Message.ViewSingleMessage;
import SQLClass.Message.ViewSingleMessageConverter;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Connection;
import java.sql.SQLException;


public class MessagesManager implements ManagerInterface{
    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;
    
    String name_class = "MessagesManager";
    
    public MessagesManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.log_manager = log_manager;
        this.connection_with_db = connection_with_db;
        this.global_var = global_var;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(SingleMessagesModel.class, new SingleMessagesModelConverter());
        builder.registerTypeAdapter(FilterMessage.class, new FilterMessageConverter());
        builder.registerTypeAdapter(ViewSingleMessage.class, new ViewSingleMessageConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        gson = builder.create();
    }

    public ResData getMessages(SingleMessagesModel item) throws MyNotAuthorizedException {
        try {
            SingleMessagesModel[] message = gson.fromJson(Select("SELECT * FROM SingleMessages " + SqlFilter(item,0,0) + ";", connection_with_db, gson), SingleMessagesModel[].class);
            for(int i = 0; i < message.length; i++) {
                MessageLinks[] message_links = gson.fromJson(Select("SELECT link_address, link_name FROM MessageLinks where user_id=" + message[i].getUser_id() + " AND message_id=" + item.getId(), connection_with_db, gson), MessageLinks[].class);
                MessageImageLinks[] image_links = gson.fromJson(Select("SELECT image_address, image_name FROM MessageImageLinks where user_id=" + message[i].getUser_id() + " AND message_id=" + item.getId(), connection_with_db, gson), MessageImageLinks[].class);
                for(int j = 0; j < message_links.length; j++) 
                    message[i].getMessage_links().add(new MessageLinks(message_links[j].getLink_address(), message_links[j].getLink_name(), message_links[j].getLink_text_link()));
                for(int j = 0; j < image_links.length; j++)
                    message[i].getImage_links().add(new MessageImageLinks(image_links[j].getImage_address(), image_links[j].getImage_name(), message_links[j].getLink_text_link()));
            }
            return new ResData(true, gson.toJson(message));
            
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getMessages(MessageModel). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getMessages(MessageModel). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getMessages(MessageModel). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getMessages(MessageModel). Exception", name_class);
        }
    }
    
    public ResData getViewMessages (FilterMessage item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT * FROM SingleMessages " + SqlFilter(item, item.getLimit(), item.getOffset()) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getViewMessages(FilterMessage). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getViewMessages(FilterMessage). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getViewMessages(FilterMessage). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getViewMessages(FilterMessage). Exception", name_class);
        }
    }
    
    public ResData getCountMessages(FilterMessage item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT COUNT(*) as result FROM SingleMessages " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getCountMessages(FilterMessage). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getCountMessages(FilterMessage). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getCountMessages(FilterMessage). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getCountMessages(FilterMessage). Exception", name_class);
        }
    }

    //Add
    
    public ResData AddMessage(SingleMessagesModel item) throws MyNotAuthorizedException {
        SQLResult res;
        String result = "", Buffressult = "";
        try {
            result = Add(item, "SingleMessages", connection_with_db, gson, true);
            res = gson.fromJson(result, SQLResult.class);
            for(int i = 0; i < item.getMessage_links().size(); i++) {
                MessageLinks pitem = new MessageLinks(res.getUser_id(),res.getId(),item.getMessage_links().get(i).getLink_address(), item.getMessage_links().get(i).getLink_name(), item.getMessage_links().get(i).getLink_text_link());
                Buffressult = Add(pitem, "MessageLinks", connection_with_db, gson, true);
            }
            for(int i = 0; i < item.getImage_links().size(); i++) {
                MessageImageLinks pitem = new MessageImageLinks(res.getUser_id(),res.getId(),item.getImage_links().get(i).getImage_address(), item.getImage_links().get(i).getImage_name(), item.getImage_links().get(i).getImage_text_link());
                Buffressult = Add(pitem, "MessageImageLinks", connection_with_db, gson, true);
            }
            return new ResData(true, result);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddMessage(MessageModel). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddMessage(MessageModel). Exeption in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddMessage(MessageModel). Exception", name_class);
        }
    }
    
    //Delete
    
    public ResData DeleteMessage(ViewSingleMessage item) throws MyNotAuthorizedException {
        try {
            //result = ExecuteQuery("DELETE FROM Messages WHERE user_id = " + item.getUser_id() + " AND invoice_id = " + item.getId() + ";", connection_with_db, gson);
            return new ResData(true, Delete(item, "SingleMessages", connection_with_db, gson));

        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteMessage(ViewMessage). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteMessage(ViewMessage). Exception", name_class);
        }
        /* catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "DeleteInvoice(ViewInvoice). Sql all delete");
        }
        **/
    }
    
    //Update
    
    public ResData UpdateMessage(SingleMessagesModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Update(item, "SingleMessages", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateMessage(MessageModel). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateMessage(MessageModel). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateMessage(MessageModel). Exception", name_class);
        }
    }
    
    public ResData Messages_select_with_count(FilterMessage item){
        try{
            ResData res, answer;
            
            answer = getViewMessages(item);
            if(!answer.getIs_success())
                return answer;
            
            res = getCountMessages(item);
            if(!res.getIs_success())
                return res;
                        
            answer.setRes("[" + answer.getRes() + "," + res.getRes().substring(1, res.getRes().length() - 1) + "]");
            answer.setIs_success(res.getIs_success());
            
            return answer;
        }
        catch(MyNotAuthorizedException ex){
            return log_manager.CallError(ex.toString(), "Messages_select_with_count(FilterMessage item). MyNotAuthorizedException", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "Messages_select_with_count(FilterMessage item). Exception", name_class);
        }
    }

    public Connection getConnection_with_db() {
        return connection_with_db;
    }

    public void setConnection_with_db(Connection connection_with_db) {
        this.connection_with_db = connection_with_db;
    }
    
    public void setErrorHttp(String ErrorHttp) {
        this.log_manager.setErrorHttp(ErrorHttp);
    }
}
