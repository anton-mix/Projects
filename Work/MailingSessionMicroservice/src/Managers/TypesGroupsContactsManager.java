package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.Contact.GroupContacts;
import SQLClass.Contact.ShortGroupContacts;
import SQLClass.Contact.ShortUserEmail;
import SQLClass.Contact.TypeGroupContacts;
import SQLClass.Contact.TypeGroupContactsConverter;
import SQLClass.GlobalVariables;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import SQLClass.Social.TelegramBots;
import SQLClass.Social.TelegramBotsConverter;
import SQLClass.Social.VkGroups;
import SQLClass.Social.VkGroupsConverter;
import SQLClass.User.UserSystemInfo;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Connection;
import java.sql.SQLException;

public class TypesGroupsContactsManager  implements ManagerInterface {
    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;
    
    String name_class = "TypesGroupsContactsManager";
    
    public TypesGroupsContactsManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.global_var = global_var;
        this.connection_with_db = connection_with_db;
        this.log_manager = log_manager;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(TypeGroupContacts.class, new TypeGroupContactsConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        builder.registerTypeAdapter(TelegramBots.class, new TelegramBotsConverter());
        builder.registerTypeAdapter(VkGroups.class, new VkGroupsConverter());
        gson = builder.create();
    }
    
    public ResData getTypesGroupsContacts (TypeGroupContacts item) throws MyNotAuthorizedException {
        try {
            TypeGroupContacts[] buff = gson.fromJson(Select("SELECT * FROM TypesGroupContacts " +  SqlFilter(item,0,0) + ";", connection_with_db, gson), TypeGroupContacts[].class);  
            
            ShortGroupContacts[] buff2 = gson.fromJson(Select("SELECT id, group_name FROM GroupsContacts WHERE user_id = '" + buff[0].getUser_id() + "' AND type_id = " + buff[0].getId() + ";", connection_with_db, gson), ShortGroupContacts[].class);
            for (int i = 0; i < buff2.length; i++) {
                buff[0].getGroups_contacts().add(buff2[i]);
            }
            
            if (item.getId() != -1) {
                ShortUserEmail[] buff3 = gson.fromJson(Select("SELECT id, email FROM UserEmails WHERE user_id = " + buff[0].getUser_id() + ";", connection_with_db, gson), ShortUserEmail[].class);
                for (int i = 0; i < buff3.length; i++) {
                    buff[0].getSender().add(buff3[i]);
                }
                
                VkGroups[] buff4 = gson.fromJson(Select("SELECT group_id, group_name, group_token FROM VkGroups WHERE user_id = " + buff[0].getUser_id() + " AND type_id = " + buff[0].getId() + ";", connection_with_db, gson), VkGroups[].class);
                if(buff4.length != 0) {
                    buff[0].setVk_group_id(buff4[0].getGroup_id());
                    buff[0].setVk_group_name(buff4[0].getGroup_name());
                    buff[0].setVk_group_token(buff4[0].getGroup_token());
                }
                
                TelegramBots[] buff5 = gson.fromJson(Select("SELECT bot_token FROM TelegramBots WHERE user_id = " + buff[0].getUser_id() + " AND type_id = " + buff[0].getId() + ";", connection_with_db, gson), TelegramBots[].class);
                if(buff5.length != 0) {
                    buff[0].setTelegram_bot_token(buff5[0].getBot_token());
                }
            }
            return new ResData(true, gson.toJson(buff));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getTypesGroupsContacts(TypeGroupContacts). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getTypesGroupsContacts(TypeGroupContacts). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getTypesGroupsContacts(TypeGroupContacts). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getTypesGroupsContacts(TypeGroupContacts). Sql select exception", name_class);
        }
    }
    
    public ResData AddTypesGroupsContacts(TypeGroupContacts item) throws MyNotAuthorizedException {
        try {
            String res = Add(item, "TypesGroupContacts", connection_with_db, gson, true);
            SQLResult res2 = gson.fromJson(res, SQLResult.class);
            String res3 = Add(new GroupContacts(item.getUser_id(),res2.getId(), item.getName(), item.getName(), 0,0,"","",item.getSender_id()), "GroupsContacts", connection_with_db, gson, true);
            
            SQLResult result = gson.fromJson(res, SQLResult.class);
            SQLResult result2 = gson.fromJson(res3, SQLResult.class);
            
            ExecuteQuery("UPDATE TypesGroupContacts SET default_group_id = " + result2.getId() + " WHERE user_id = " + item.getUser_id() + " AND id = " + result.getId() + ";", connection_with_db, gson);
            
            for(int i = 0; i < item.getGroups_contacts().size();i++) {
                ExecuteQuery("UPDATE GroupsContacts SET type_id = " + result.getId() + " WHERE user_id = " + item.getUser_id() + " AND id = " + item.getGroups_contacts().get(i).getId() + ";", connection_with_db, gson);
            }
            if(!"".equals(item.getVk_group_token())) {
                String res4 = Add(new VkGroups(item.getUser_id(), res2.getId(), item.getVk_group_name(), item.getVk_group_token(), item.getVk_group_id(), 1,0), "VkGroups", connection_with_db, gson, true);
            }
            
            if(!"".equals(item.getTelegram_bot_token())) {
                String res5 = Add(new TelegramBots(item.getUser_id(), res2.getId(), item.getTelegram_bot_token(), "", 1, 0), "TelegramBots", connection_with_db, gson, true);
            }
            
            return new ResData(true, res);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddTypesGroupsContacts(TypeGroupContacts). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddTypesGroupsContacts(TypeGroupContacts). Exception in parse", "TypesGroupsContactsManagers");
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddTypesGroupsContacts(TypeGroupContacts). Exception", "TypesGroupsContactsManagers");
        }
    }
    
    public ResData DeleteTypesGroupsContacts(TypeGroupContacts item) throws MyNotAuthorizedException {
        try {   
            UserSystemInfo buff2 = gson.fromJson(getUserSystemInfo(gson.fromJson("{\"user_id\":\"" + item.getUser_id() + "\"}", UserSystemInfo.class)).getRes(), UserSystemInfo[].class)[0];
            if (buff2.getType_group_contacts_id() != item.getId()) {
                ExecuteQuery("Delete from GroupsContacts WHERE user_id = " + item.getUser_id() + " AND type_id = " + item.getId() + ";", connection_with_db, gson);
                return new ResData(true, Delete(item, "TypesGroupContacts", connection_with_db, gson));
            }
            else {
                return new ResData(true, gson.toJson(new SQLResult("Ok"), SQLResult.class));
            }
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteTypesGroupsContacts(TypeGroupContacts). Sql delete", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "DeleteTypesGroupsContacts(TypeGroupContacts). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteTypesGroupsContacts(TypeGroupContacts). exception", name_class);
        }
    }
    
    public ResData UpdateTypesGroupsContacts(TypeGroupContacts item) throws MyNotAuthorizedException {
        try {
            String res = "",res2 = "";
            res = Update(item, "TypesGroupContacts", connection_with_db, gson);
            for(int i = 0; i < item.getGroups_contacts().size(); i++) {
                res2 = getTypesGroupsContacts(gson.fromJson("{\"user_id\": " + item.getUser_id() + ", \"default_group_id\": " + item.getGroups_contacts().get(i).getId() + "}", TypeGroupContacts.class)).getRes();
                if("[]".equals(res2))
                    ExecuteQuery("UPDATE GroupsContacts SET type_id = " + item.getId() + " WHERE user_id = '" + item.getUser_id() + "' AND id = " + item.getGroups_contacts().get(i).getId() + ";", connection_with_db, gson);
            }         
            VkGroups[] vkGroups = null;
            try {
                vkGroups = gson.fromJson(Select("SELECT * FROM VkGroups  WHERE user_id = " + item.getUser_id() + " AND type_id = '" + item.getId() + "';", connection_with_db, gson), VkGroups[].class);
            } catch (MyTypeException ex) {
                return log_manager.CallError(ex.toString(), "UpdateTypesGroupsContacts(TypeGroupContacts). getVGgroup", name_class);
            }
            
            String res4 = "";
            if(!"".equals(item.getVk_group_token())) {
                if(vkGroups != null && vkGroups.length > 0) {
                    res4 = ExecuteQuery("UPDATE VkGroups SET group_name  = '" + item.getVk_group_name() + "', group_token = '" + item.getVk_group_token() + "', group_id  = " + item.getVk_group_id() + " WHERE user_id = '" + item.getUser_id() + "' AND type_id = " + item.getId() + ";", connection_with_db, gson);
                }
                else  {
                    res4 = Add(new VkGroups(item.getUser_id(), item.getId(), item.getVk_group_name(), item.getVk_group_token(), item.getVk_group_id(), 1,0), "VkGroups", connection_with_db, gson, true);   
                }
            }
            else {
                if(vkGroups != null && vkGroups.length > 0) {
                    res4 = ExecuteQuery("DELETE FROM VkGroups WHERE user_id = '" + item.getUser_id() + "' AND type_id = " + item.getId() + ";", connection_with_db, gson);
                }
            }
            
            TelegramBots[] telegramBots = null;
            try {
                telegramBots = gson.fromJson(Select("SELECT * FROM TelegramBots WHERE user_id = " + item.getUser_id() + " AND type_id = '" + item.getId() + "';", connection_with_db, gson), TelegramBots[].class);
            } catch (MyTypeException ex) {
                return log_manager.CallError(ex.toString(), "UpdateTypesGroupsContacts(TypeGroupContacts). getVGgroup", name_class);
            }
            
            String res5 = "";
            if(!"".equals(item.getTelegram_bot_token())) {
                if(telegramBots != null && telegramBots.length > 0) { 
                    res5 = ExecuteQuery("UPDATE TelegramBots SET bot_token = '" + item.getTelegram_bot_token()+ "' WHERE user_id = '" + item.getUser_id() + "' AND type_id = " + item.getId() + ";", connection_with_db, gson);
                }
                else {
                    res5 = Add(new TelegramBots(item.getUser_id(), item.getId(), item.getTelegram_bot_token(), "", 1, 0), "TelegramBots", connection_with_db, gson, true);
                }
            }
            else 
                if(telegramBots != null && telegramBots.length > 0) {
                    res5 = ExecuteQuery("DELETE FROM TelegramBots WHERE user_id = '" + item.getUser_id() + "' AND type_id = " + item.getId() + ";", connection_with_db, gson);
                }
                
            return new ResData(true, res);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateTypesGroupsContacts(TypeGroupContacts). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateTypesGroupsContacts(TypeGroupContacts). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateTypesGroupsContacts(TypeGroupContacts). Exception", name_class);
        }
    }
    
    public ResData getUserSystemInfo(UserSystemInfo item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT * FROM UserSystemInfo " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getUserSystemInfo(UserSystemInfo). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getUserSystemInfo(UserSystemInfo). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getUserSystemInfo(UserSystemInfo). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getUserSystemInfo(UserSystemInfo). Exception", name_class);
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
