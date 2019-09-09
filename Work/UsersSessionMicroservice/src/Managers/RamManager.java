package Managers;

import MyEnum.RestRamEnum;
import QueueClass.ArrUserPasswordKeyValue;
import QueueClass.ResQueueItem;
import QueueClass.UserAuthModel;
import SQLClass.SQLResult;
import com.google.gson.Gson;
import java.io.IOException;

public class RamManager {
    LogManager log_manager;
    MessageManager message_manager;
    Gson gson;
    
    String name_class = "RamManager";
    
    public RamManager(LogManager log_manager, MessageManager message_manager, Gson gson) {
        this.log_manager = log_manager;
        this.message_manager = message_manager;
        this.gson = gson;
    }
    
    public SQLResult CheckSessionKey(String SessionKey) {
        String res = "";
        try {
            res = message_manager.call("{}", gson, SessionKey, RestRamEnum.Check_session_key.ordinal(), "ram_queue").getRes();
        } catch (IOException ex) {
            res = log_manager.CallError(ex.toString(), "CheckSessionKey(String). IOException", name_class).getRes();
        } catch (InterruptedException ex) {
            res = log_manager.CallError(ex.toString(), "CheckSessionKey(String). InterruptedException", name_class).getRes();
        } catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "CheckSessionKey(String). Exeption", name_class).getRes();
        } 
        return gson.fromJson(res, SQLResult.class);
    }
    
    public SQLResult AddSessionKey(int id) {
        String res = "";
        try {
            res = message_manager.call(new ResQueueItem(id, ""), gson, "", RestRamEnum.Add_session_key.ordinal(), "ram_queue").getRes();
        } catch (IOException ex) {
            res = log_manager.CallError(ex.toString(), "AddSessionKey(int). IOException", name_class).getRes();
        } catch (InterruptedException ex) {
            res = log_manager.CallError(ex.toString(), "AddSessionKey(int). InterruptedException", name_class).getRes();
        } catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "AddSessionKey(int). Exeption", name_class).getRes();
        } 
        return gson.fromJson(res, SQLResult.class);
    }
    
    public SQLResult CheckAndDeleteUserAuthMessage(String SessionKey) {
        String res = "";
        try {
            res = message_manager.call("{}", gson, SessionKey, RestRamEnum.Check_and_delete_user_auth_message.ordinal(), "ram_queue").getRes();
        } catch (IOException ex) {
            res = log_manager.CallError(ex.toString(), "CheckAndDeleteUserAuthMessage(String). IOException", name_class).getRes();
        } catch (InterruptedException ex) {
            res = log_manager.CallError(ex.toString(), "CheckAndDeleteUserAuthMessage(String). InterruptedException", name_class).getRes();
        } catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "CheckAndDeleteUserAuthMessage(String). Exeption", name_class).getRes();
        } 
        return gson.fromJson(res, SQLResult.class);
    }
            
    public SQLResult CheckAuthMessageKeyEmail(String Email) {
        String res = "";
        try {
            res = message_manager.call("{}", gson, Email, RestRamEnum.Check_auth_message_key_email.ordinal(), "ram_queue").getRes();
        } catch (IOException ex) {
            res = log_manager.CallError(ex.toString(), "CheckAuthMessageKeyEmail(String). Exeption", name_class).getRes();
        } catch (InterruptedException ex) {
            res = log_manager.CallError(ex.toString(), "CheckAuthMessageKeyEmail(String). Exeption", name_class).getRes();
        } catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "CheckAuthMessageKeyEmail(String). Exeption", name_class).getRes();
        } 
        return gson.fromJson(res, SQLResult.class);
    }
            
    public SQLResult Add_auth_message_key(UserAuthModel item) {
        String res = "";
        try {
            res = message_manager.call(gson.toJson(item, UserAuthModel.class), gson, "", RestRamEnum.Add_auth_message_key.ordinal(), "ram_queue").getRes();
        } catch (IOException ex) {
            res = log_manager.CallError(ex.toString(), "Add_auth_message_key(String). Exeption", name_class).getRes();
        } catch (InterruptedException ex) {
            res = log_manager.CallError(ex.toString(), "Add_auth_message_key(String). Exeption", name_class).getRes();
        } catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "Add_auth_message_key(String). Exeption", name_class).getRes();
        } 
        return gson.fromJson(res, SQLResult.class);
    }

    public SQLResult CheckChangePassMessageKeyLogin(String Login) {
        String res = "";
        try {
            res = message_manager.call("{}", gson, Login, RestRamEnum.Check_change_pass_message_key_login.ordinal(), "ram_queue").getRes();
        } catch (IOException ex) {
            res = log_manager.CallError(ex.toString(), "CheckChangePassMessageKeyLogin(String). Exeption", name_class).getRes();
        } catch (InterruptedException ex) {
            res = log_manager.CallError(ex.toString(), "CheckChangePassMessageKeyLogin(String). Exeption", name_class).getRes();
        } catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "CheckChangePassMessageKeyLogin(String). Exeption", name_class).getRes();
        }
        return gson.fromJson(res, SQLResult.class);
    }
    
    public SQLResult AddChangePassMessageKey(String login) {
        String res = "";
        try {
            res = message_manager.call(gson.toJson(new ArrUserPasswordKeyValue(login), ArrUserPasswordKeyValue.class), gson, "", RestRamEnum.Add_change_pass_message_key.ordinal(), "ram_queue").getRes();
        } catch (IOException ex) {
            res = log_manager.CallError(ex.toString(), "AddChangePassMessageKey(String). Exeption", name_class).getRes();
        } catch (InterruptedException ex) {
            res = log_manager.CallError(ex.toString(), "AddChangePassMessageKey(String). Exeption", name_class).getRes();
        } catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "AddChangePassMessageKey(String). Exeption", name_class).getRes();
        } 
        return gson.fromJson(res, SQLResult.class);
    }
    
    public SQLResult CheckAndDeleteChangePassMessageKey(String SessionKey) {
        String res = "";
        try {
            res = message_manager.call("{}", gson, SessionKey, RestRamEnum.Check_and_delete_change_pass_message_key.ordinal(), "ram_queue").getRes();
            System.out.println(res);
        } catch (IOException ex) {
            res = log_manager.CallError(ex.toString(), "CheckAndDeleteChangePassMessageKey(String). Exeption", name_class).getRes();
        } catch (InterruptedException ex) {
            res = log_manager.CallError(ex.toString(), "CheckAndDeleteChangePassMessageKey(String). Exeption", name_class).getRes();
        } catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "CheckAndDeleteChangePassMessageKey(String). Exeption", name_class).getRes();
        } 
        return gson.fromJson(res, SQLResult.class);
    }
    
    public SQLResult DeleteChangePassMessageItemFromLogin(String login) {
        String res = "";
        try {
            res = message_manager.call(gson.toJson(new ArrUserPasswordKeyValue(login), ArrUserPasswordKeyValue.class), gson, "", RestRamEnum.Delete_change_pass_message_item_from_login.ordinal(), "ram_queue").getRes();
        } catch (IOException ex) {
            res = log_manager.CallError(ex.toString(), "CheckAndDeleteChangePassMessageKey(String). Exeption", name_class).getRes();
        } catch (InterruptedException ex) {
            res = log_manager.CallError(ex.toString(), "CheckAndDeleteChangePassMessageKey(String). Exeption", name_class).getRes();
        } catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "CheckAndDeleteChangePassMessageKey(String). Exeption", name_class).getRes();
        } 
        return gson.fromJson(res, SQLResult.class);
    }
    
    public SQLResult DeleteUserAuthMessageItemForEmail(String email) {
        String res = "";
        try {
            res = message_manager.call(gson.toJson(new UserAuthModel(email), UserAuthModel.class), gson, "", RestRamEnum.Delete_user_auth_message_item_for_email.ordinal(), "ram_queue").getRes();
        } catch (IOException ex) {
            res = log_manager.CallError(ex.toString(), "DeleteUserAuthMessageItemForEmail(String). Exeption", name_class).getRes();
        } catch (InterruptedException ex) {
            res = log_manager.CallError(ex.toString(), "DeleteUserAuthMessageItemForEmail(String). Exeption", name_class).getRes();
        } catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "DeleteUserAuthMessageItemForEmail(String). Exeption", name_class).getRes();
        } 
        return gson.fromJson(res, SQLResult.class);
    }
}