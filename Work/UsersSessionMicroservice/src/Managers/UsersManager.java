package Managers;

import MessageClass.MessageAuth;
import MessageClass.MessageChangePassword;
import MessageClass.MessageChangePasswordConverter;
import MyEnum.MailMessageEnum;
import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import QueueClass.ArrUserPasswordKeyValue;
import QueueClass.ArrUserPasswordKeyValueConverter;
import QueueClass.ResQueueItem;
import QueueClass.UserAuthModel;
import QueueClass.UserUpdatePassword;
import RestObject.RessiveClass;
import RestObject.RessiveClassConverter;
import SQLClass.Contact.GroupContacts;
import SQLClass.Contact.TypeGroupContacts;
import SQLClass.GlobalVariables;
import SQLClass.HttpSettings.UserEmails;
import SQLClass.Product.GroupProducts;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import SQLClass.User.UserModel;
import SQLClass.User.UserModelConverter;
import SQLClass.User.UserSystemInfo;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import java.io.IOException;

import java.sql.Connection;
import java.sql.SQLException;
import java.util.logging.Level;
import java.util.logging.Logger;

public class UsersManager implements ManagerInterface  {
    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;
    JsonParser parser = new JsonParser();
    String name_class = "UsersManager";
    
    public UsersManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.log_manager = log_manager;
        this.connection_with_db = connection_with_db;
        this.global_var = global_var;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(UserModel.class, new UserModelConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        builder.registerTypeAdapter(MessageChangePassword.class, new MessageChangePasswordConverter());
        builder.registerTypeAdapter(RessiveClass.class, new RessiveClassConverter());
        builder.registerTypeAdapter(ArrUserPasswordKeyValue.class, new ArrUserPasswordKeyValueConverter());
        gson = builder.create();
    }

    //Select
    
    public ResData getUsers() {
        try {
            return new ResData(true, Select("SELECT * FROM Users ;", connection_with_db, gson));   
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getUsers(). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getUsers(). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getUsers(). Exception", name_class);
        }
    }
    
 
    
    public ResData getUser (UserModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT * FROM Users " + SqlFilter(item,0,0) + ";", connection_with_db, gson));   
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getUser(User). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getUser(User). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getUser(User). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getUser(User). Exception", name_class);
        }
    }
    
    public ResData getOneUser(UserModel item, RamManager ram_manager) throws MyNotAuthorizedException {
        try {
            String res = "";
            if ("".equals(item.getUsername()) || "".equals(item.getPassword())) {
                return log_manager.CallError("Field username or password is null", "call to queue", name_class);
            } else {
                res = Select("SELECT * FROM Users WHERE username = '" + item.getUsername() + "' AND password = '" + item.getPassword() + "';", connection_with_db, gson);
            }

            if ("[]".equals(res)) {
                return new ResData(true, gson.toJson(new SQLResult("User not found")));
            }

            if (res.charAt(0) == '[') {
                JsonObject json_obj = parser.parse(res).getAsJsonArray().get(0).getAsJsonObject();
                SQLResult resQueue = ram_manager.AddSessionKey(json_obj.get("id").getAsInt());
                if (resQueue.getId() == -1) {
                    return new ResData(false, gson.toJson(resQueue, SQLResult.class));
                }
                ResQueueItem buff = gson.fromJson(resQueue.getResult(), ResQueueItem.class);
                json_obj.addProperty("user_id", buff.getSession_key());
                return new ResData(true, "[" + gson.toJson(json_obj) + "]");
            } else {
                return new ResData(true, res);
            } 
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getOneUser(UserModel). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getOneUser(UserModel). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getOneUser(UserModel). Exception", name_class);
        }
    }
    
    public ResData getUserForLogin (String username) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT * FROM Users WHERE username = '" + username + "';", connection_with_db, gson));   
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getUserForLogin(String). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getUserForLogin(String). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getUserForLogin(String). Exception", name_class);
        }
    }
    
    public ResData getUserEmails(UserEmails item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT * FROM UserEmails WHERE user_id = " + item.getId() + ";", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getUserEmails(UserModel). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getUserEmails(UserModel). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getUserEmails(UserModel). Exception", name_class);
        }
    }
    
    public ResData getUserRegEmail(UserModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT id, email FROM Users WHERE id = " + item.getId() + ";", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getUserRegEmail(UserModel). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getUserRegEmail(UserModel). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getUserRegEmail(UserModel). Exception", name_class);
        }
    }
    
    //Add
    
    public ResData AddUser(UserAuthModel item, RamManager ram_manager, UserSystemInfoManager userSystemInfo_manager) throws MyNotAuthorizedException {
        try {
            SQLResult resQueue, res;
            resQueue = ram_manager.CheckAndDeleteUserAuthMessage(item.getSession_id());
            if (resQueue.getId() == -1) {
                return new ResData(false, gson.toJson(resQueue, SQLResult.class));
            }
            UserAuthModel buff2 = gson.fromJson(resQueue.getResult(), UserAuthModel.class);

            if ("".equals(buff2.getEmail())) {
                return new ResData(false, gson.toJson(new SQLResult("Key not found")));
            } else {
                item.AddFilds(buff2);
                UserModel buff3 = new UserModel("", "", "", "", "", "", "", "");
                buff3.addFilds(item);

                res = gson.fromJson(Add(buff3, "Users", connection_with_db, gson, false), SQLResult.class);
                resQueue = gson.fromJson(Execute("call create_default_user_data(" + res.getId() + ")", connection_with_db, gson), SQLResult.class);
                
                resQueue = ram_manager.AddSessionKey(res.getId());
                if (resQueue.getId() == -1) {
                    return new ResData(false, gson.toJson(resQueue, SQLResult.class));
                }
                JsonObject json_obj = parser.parse(gson.toJson(res, SQLResult.class)).getAsJsonObject();
                json_obj.addProperty("user_id", resQueue.getResult());
                return new ResData(true, gson.toJson(json_obj)); 
            }
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddUser(UserModel). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddUser(UserModel). Exeption in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddUser(UserModel). Exception", name_class);
        }
    }
    
    //Delete
    
    public ResData DeleteUser(UserModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Delete(item, "Users", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteUser(UserModel). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteUser(UserModel). Exception", name_class);
        }
    }
    
    //Update
    
    public ResData UpdateUser(UserModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Update(item, "Users", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateUser(UserModel). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateUser(UserModel). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateUser(UserModel). Exception", name_class);
        }
    }
    
    //Other function
    
    public ResData UserCancel(UserAuthModel item, RamManager ram_manager) {
        SQLResult resQueue;
        try {
            resQueue = ram_manager.CheckAndDeleteUserAuthMessage(item.getSession_id());
            if (resQueue.getId() == -1) {
                return new ResData(false, gson.toJson(resQueue, SQLResult.class));
            }

            UserAuthModel buff2 = gson.fromJson(resQueue.getResult(), UserAuthModel.class);

            if ("".equals(buff2.getEmail())) {
                return new ResData(false, gson.toJson(new SQLResult("Key not found")));
            } else {
                return new ResData(true, gson.toJson(new SQLResult("Cancel success")));
            }
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UserCancel(UserAuthModel, RamManager). Exception", name_class);
        }
    }
    
    public ResData UserPreReg(UserAuthModel item, RamManager ram_manager, MessageManager message_manager) {
        SQLResult resQueue;
        try {
            resQueue = ram_manager.CheckAuthMessageKeyEmail(item.getEmail());
            if (resQueue.getId() == -1) {
                return new ResData(false, gson.toJson(resQueue, SQLResult.class));
            }
            
            if (resQueue.getId() == 1) {

                if ("not Ok".equals(resQueue.getResult())) {
                    
                    resQueue = ram_manager.DeleteUserAuthMessageItemForEmail(item.getEmail());
                    if (resQueue.getId() == -1) {
                        return new ResData(false, gson.toJson(resQueue, SQLResult.class));
                    }
                }
                resQueue = ram_manager.Add_auth_message_key(item);
                if (resQueue.getId() == -1) {
                    return new ResData(false, gson.toJson(resQueue, SQLResult.class));
                }
                JsonObject json_obj = parser.parse(gson.toJson(new MessageAuth(item.getName(), resQueue.getResult(), item.getEmail(), "info@utp1.ru"), MessageAuth.class)).getAsJsonObject();
                return message_manager.call(gson.toJson(new RessiveClass(json_obj, MailMessageEnum.Registr_user.ordinal()), RessiveClass.class), "mail_queue");
            }
            else {
                return new ResData(true, gson.toJson(new SQLResult("Ok"), SQLResult.class));
            }
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UserPreReg(UserAuthModel, RamManager, MessageManager). Exception", name_class);
        }
    }
    
    public ResData UserCheckParam(UserModel item) {
        try {
            ResData response = getUser(item);

            if (response.getIs_success() == false) {
                return response;
            }

            if (response.getRes().charAt(0) != '[') {
                response = new ResData(true, gson.toJson(new SQLResult("ne Ok")));
            } else {
                if ("[]".equals(response.getRes())) {
                    response = new ResData(true, gson.toJson(new SQLResult("Ok")));
                } else {
                    response = new ResData(false, gson.toJson(new SQLResult("ne Ok")));
                }
            }
            return response;
        }
        catch (MyNotAuthorizedException ex) {
            return log_manager.CallError(ex.toString(), "UserCheckParam(UserModel). MyNotAuthorizedException", name_class);
        }
        catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UserCheckParam(UserModel). Exception", name_class);
        }
    }
    
    public ResData UserBeginChangePassword(UserModel item, RamManager ram_manager, MessageManager message_manager) {
        SQLResult resQueue;
        ResData response;
        try {
            response = getUserForLogin(item.getUsername());
            if (response.getIs_success() == false) {
                return response;
            }

            if ("[]".equals(response.getRes())) {
                return new ResData(false, gson.toJson(new SQLResult("User not found")));
            } else if (response.getRes().charAt(0) == '[') {
                UserModel buff = gson.fromJson(response.getRes(), UserModel[].class)[0];

                resQueue = ram_manager.CheckChangePassMessageKeyLogin(buff.getUsername());
                if (resQueue.getId() == -1) {
                    return new ResData(false, gson.toJson(resQueue, SQLResult.class));
                }
                if (resQueue.getId() == 1) {
                    if ("not Ok".equals(resQueue.getResult())) {
                        resQueue = ram_manager.DeleteChangePassMessageItemFromLogin(buff.getUsername());
                        if (resQueue.getId() == -1) {
                            return new ResData(false, gson.toJson(resQueue, SQLResult.class));
                        }
                    }

                    resQueue = ram_manager.AddChangePassMessageKey(buff.getUsername());
                    if (resQueue.getId() == -1) {
                        return new ResData(false, gson.toJson(resQueue, SQLResult.class));
                    }
                    JsonObject json_obj = parser.parse(gson.toJson(new MessageChangePassword(resQueue.getResult(), buff.getEmail(), "info@utp1.ru"), MessageChangePassword.class)).getAsJsonObject();
                    response = message_manager.call(gson.toJson(new RessiveClass(json_obj, MailMessageEnum.SendChangePassword.ordinal()), RessiveClass.class), "mail_queue");

                } else {
                    response = new ResData(false, gson.toJson(new SQLResult("Error in ram_microservice")));
                }
            }
            return response;
        } catch (IOException ex) {
            return log_manager.CallError(ex.toString(), "UserCheckParam(UserModel). IOException", name_class);
        } catch (MyNotAuthorizedException ex) {
            return log_manager.CallError(ex.toString(), "UserCheckParam(UserModel). MyNotAuthorizedException", name_class);
        } catch (InterruptedException ex) {
            return log_manager.CallError(ex.toString(), "UserCheckParam(UserModel). InterruptedException", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UserCheckParam(UserModel). Exception", name_class);
        }
    }
    
    public ResData UserEndChangePassword(UserModel item, RamManager ram_manager) {
        ResData response = new ResData();
        SQLResult resQueue;
        try {

            resQueue = ram_manager.CheckAndDeleteChangePassMessageKey(item.getUsername());
            if (resQueue.getId() == -1) {
                return new ResData(false, gson.toJson(resQueue, SQLResult.class));
            }
        
            if ("not Ok".equals(resQueue.getResult())) {
                response = new ResData(true, gson.toJson(new SQLResult("Key not found")));
            } else {
                ArrUserPasswordKeyValue buff2 = gson.fromJson(resQueue.getResult(), ArrUserPasswordKeyValue.class);
                response = new ResData(true, ExecuteQuery("UPDATE Users set password = '" + item.getPassword() + "' where username = '" + buff2.getLogin() + "';", connection_with_db, gson));
            }
            return response;
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UserEndChangePassword(UserUpdatePassword). SQLException", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UserEndChangePassword(UserUpdatePassword). MyParseException", name_class);
        } catch (MyNotAuthorizedException ex) {
            return log_manager.CallError(ex.toString(), "UserEndChangePassword(UserUpdatePassword). MyNotAuthorizedException", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UserEndChangePassword(UserUpdatePassword). Exception", name_class);
        }
    }
    
    public ResData UserUpdateAccountPassword(UserUpdatePassword item) {
        ResData response;
        try {
            response = new ResData(true, Select("select * from Users where id = " + item.getUser_id() + " and password = '" + item.getOld_password() + "';", connection_with_db, gson));
            if (response.getIs_success() == false) {
                return response;
            }
            if ("[]".equals(response.getRes())) {
                response = new ResData(false, gson.toJson(new SQLResult("Invalid password")));
            } else if (response.getRes().charAt(0) == '[') {
                response = new ResData(true, ExecuteQuery("UPDATE Users set password = '" + item.getNew_password() + "' where id = " + item.getUser_id() + ";", connection_with_db, gson));
            }
            return response;
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UserEndChangePassword(UserUpdatePassword). SQLException", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UserEndChangePassword(UserUpdatePassword). MyParseException", name_class);
        } catch (MyNotAuthorizedException ex) {
            return log_manager.CallError(ex.toString(), "UserEndChangePassword(UserUpdatePassword). MyNotAuthorizedException", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "UserEndChangePassword(UserUpdatePassword). MyTypeException", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UserEndChangePassword(UserUpdatePassword). Exception", name_class);
        }
    }
    
    public ResData AddGroupsProducts(GroupProducts item) throws MyNotAuthorizedException {
        try {
            String res = Add(item, "GroupsProducts", connection_with_db, gson, true);
            int id = gson.fromJson(res, SQLResult.class).getId();
            for(int i = 0; i < item.getProducts().size();i++) {
                res = ExecuteQuery("UPDATE Products SET group_id = " + id + " WHERE user_id = " + item.getUser_id() + " AND id = " + item.getProducts().get(i).getId() + ";", connection_with_db, gson);
            }
            return new ResData(true, res);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddGroupsProducts(GroupProducts). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddGroupsProducts(GroupProducts). Exeption in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddGroupsProducts(GroupProducts). Exception", name_class);
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
            
            return new ResData(true, res);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddTypesGroupsContacts(TypeGroupContacts). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddTypesGroupsContacts(TypeGroupContacts). Exeption in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddTypesGroupsContacts(TypeGroupContacts). Exception", name_class);
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
