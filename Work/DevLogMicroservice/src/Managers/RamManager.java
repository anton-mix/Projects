package Managers;

import MyEnum.RestRamEnum;
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
            res = message_manager.call("{}", gson, SessionKey, RestRamEnum.Check_session_key.ordinal(), RestRamEnum.Check_session_key.name(), "ram_queue");
        } catch (IOException ex) {
            res = log_manager.CallError(ex.toString(), "CheckSessionKey(String). IOException", name_class);
        } catch (InterruptedException ex) {
            res = log_manager.CallError(ex.toString(), "CheckSessionKey(String). InterruptedException", name_class);
        } catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "CheckSessionKey(String). Exeption", name_class);
        } 
        return gson.fromJson(res, SQLResult.class);
    }
}