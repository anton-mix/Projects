package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.Contact.ContactsToGroupContact;
import SQLClass.Contact.GroupContacts;
import SQLClass.Contact.GroupContactsConverter;
import SQLClass.Contact.TypeGroupContacts;
import SQLClass.Contact.ViewGroupContacts;
import SQLClass.Contact.ViewGroupContactsConverter;
import SQLClass.Contact.ViewTypeGroupContacts;
import SQLClass.Contact.ViewTypeGroupContactsConverter;
import SQLClass.GlobalVariables;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.sql.Connection;
import java.sql.SQLException;

public class GroupsContactsManager implements ManagerInterface {

    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;

    String name_class = "GroupsContactsManager";

    public GroupsContactsManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.log_manager = log_manager;
        this.connection_with_db = connection_with_db;
        this.global_var = global_var;

        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(ViewGroupContacts.class, new ViewGroupContactsConverter());
        builder.registerTypeAdapter(GroupContacts.class, new GroupContactsConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        builder.registerTypeAdapter(ViewTypeGroupContacts.class, new ViewTypeGroupContactsConverter());
        gson = builder.create();
    }

    public ResData getGroupContacts(ViewGroupContacts item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT * FROM GroupsContacts " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getGroupsContact(GroupContactsView). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getGroupsContact(GroupContactsView). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getGroupsContact(GroupContactsView). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getGroupsContact(GroupContactsView). Exception", name_class);
        }
    }

    public ResData getGroupsContactsView(ViewGroupContacts item) throws MyNotAuthorizedException {
        try {
            ViewGroupContacts[] buff = gson.fromJson(Select("SELECT * FROM GroupsContacts " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson), ViewGroupContacts[].class);
            for (int i = 0; i < buff.length; i++) {
                TypeGroupContacts[] buff2 = gson.fromJson(Select("SELECT name FROM TypesGroupContacts  WHERE user_id = " + buff[i].getUser_id() + " AND id = " + buff[i].getType_id() + ";", connection_with_db, gson), TypeGroupContacts[].class);
                buff[i].setType_name(buff2[0].getName());
            }
            return new ResData(true, gson.toJson(buff));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getGroupsContactsView(GroupContactsView). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getGroupsContactsView(GroupContactsView). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getGroupsContactsView(GroupContactsView). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getGroupsContactsView(GroupContactsView). Exception", name_class);
        }
    }

    public ResData AddGroupsContacts(GroupContacts item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Add(item, "GroupsContacts", connection_with_db, gson, true));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddGroupsContacts(GroupContacts). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddGroupsContacts(GroupContacts). Exception in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddGroupsContacts(GroupContacts). Exception", name_class);
        }
    }

    public ResData DeleteGroupsContacts(GroupContacts item) throws MyNotAuthorizedException {
        try {
            String res = Select("SELECT * FROM TypesGroupContacts WHERE default_group_id = " + item.getId() + " AND user_id = " + item.getUser_id() + ";", connection_with_db, gson);
            if ("[]".equals(res)) {
                res = Delete(item, "GroupsContacts", connection_with_db, gson);
            } else {
                res = gson.toJson(new SQLResult("Ok"), SQLResult.class);
            }
            return new ResData(true, res);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteGroupsContacts(GroupContacts). Sql delete", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "DeleteGroupsContacts(GroupContacts). Sql delete my type", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteGroupsContacts(GroupContacts). Sql delete exception", name_class);
        }
    }

    public ResData UpdateGroupsContacts(GroupContacts item) throws MyNotAuthorizedException {
        try {
            String res = Select("SELECT * FROM TypesGroupContacts WHERE default_group_id = " + item.getId() + " AND user_id = " + item.getUser_id() + ";", connection_with_db, gson);
            if ("[]".equals(res)) {
                res = Update(item, "GroupsContacts", connection_with_db, gson);
            } else {
                res = gson.toJson(new SQLResult("Ok"), SQLResult.class);
            }
            return new ResData(true, res);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateGroupsProducts(GroupProducts). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateGroupsProducts(GroupProducts). Exception In parse", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "UpdateGroupsProducts(GroupProducts). Exception In Type", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateGroupsProducts(GroupProducts). Exception", name_class);
        }
    }

    public ResData AddContactToGroupContacts(ContactsToGroupContact item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Add(item, "Contacts_to_GroupContacts", connection_with_db, gson, true));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddContactToGroupContacts(ContactsToGroupContact). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddContactToGroupContacts(ContactsToGroupContact). Exception in parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddContactToGroupContacts(ContactsToGroupContact). Exception", name_class);
        }
    }

    public ResData DeleteContactToGroupContacts(ContactsToGroupContact item) throws MyNotAuthorizedException {
        try {
            return new ResData(true,  Delete(item, "Contacts_to_GroupContacts", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteContactToGroupContacts. Sql add", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteContactToGroupContacts. Sql add exception", name_class);
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
