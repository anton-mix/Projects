package Managers;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import MyInterface.ManagerInterface;
import SQLClass.Contact.ContactEmail;
import SQLClass.Contact.ContactModel;
import SQLClass.Contact.ContactModelConverter;
import SQLClass.Contact.ContactPhone;
import SQLClass.Contact.ContactSocial;
import SQLClass.Contact.ContactsToGroupContact;
import SQLClass.Contact.FilterContact;
import SQLClass.Contact.GroupContacts;
import SQLClass.Contact.GroupContactsConverter;
import SQLClass.Contact.ShortGroupContacts;
import SQLClass.Contact.TypeGroupContacts;
import SQLClass.Contact.ViewContact;
import SQLClass.Contact.ViewContactConverter;
import SQLClass.GlobalVariables;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import SQLClass.User.UserSystemInfo;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import java.sql.Connection;
import java.sql.SQLException;
import java.util.ArrayList;

public class ContactsManager implements ManagerInterface {

    LogManager log_manager;
    GlobalVariables global_var;
    Connection connection_with_db = null;
    Gson gson;
    String name_class = "ContactsManager";

    public ContactsManager(GlobalVariables global_var, Connection connection_with_db, LogManager log_manager) {
        this.log_manager = log_manager;
        this.connection_with_db = connection_with_db;
        this.global_var = global_var;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(ContactModel.class, new ContactModelConverter());
        builder.registerTypeAdapter(GroupContacts.class, new GroupContactsConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        builder.registerTypeAdapter(ViewContact.class, new ViewContactConverter());
        gson = builder.create();
    }

    public ResData getContact(ViewContact item) throws MyNotAuthorizedException {
        try {
            ContactModel[] contact = gson.fromJson(Select("SELECT * FROM Contacts " + SqlFilter(item,0,0) + ";", connection_with_db, gson), ContactModel[].class);
            
            ContactEmail[] buff2 = gson.fromJson(Select("SELECT * FROM ContactEmails WHERE user_id = " + item.getUser_id() + " AND contact_id = " + contact[0].getId()+ ";", connection_with_db, gson),ContactEmail[].class);
            for(int i = 0; i < buff2.length; i++) {
                contact[0].getEmails().add(buff2[i].getEmail());
            }
            
            ContactPhone[] buff3 = gson.fromJson(Select("SELECT * FROM ContactPhones WHERE user_id = " + item.getUser_id() + " AND contact_id = " + contact[0].getId()+ ";", connection_with_db, gson),ContactPhone[].class);
            for(int i = 0; i < buff3.length; i++) {
                contact[0].getPhones().add(buff3[i].getPhone());
            }
            
            ContactSocial[] buff4 = gson.fromJson(Select("SELECT * FROM ContactSocial WHERE user_id = '" + item.getUser_id() + "' AND contact_id = " + contact[0].getId()+ ";", connection_with_db, gson), ContactSocial[].class);
            for(int i = 0; i < buff4.length; i++) {
                contact[0].getSocial().add(buff4[i].getSocial());
            }
            
            GroupContacts[] buff5 = gson.fromJson(Select("SELECT c.id, c.group_name FROM Contacts a, Contacts_to_GroupContacts b, GroupsContacts c WHERE c.user_id = '" + item.getUser_id() + "' AND c.id = b.group_contacts_id AND a.id = b.contacts_id  and a.id = " + item.getId() + ";", connection_with_db, gson), GroupContacts[].class);
            for(int i = 0; i < buff5.length; i++) {
                contact[0].getGroups_contacts().add(new ShortGroupContacts(buff5[i].getId(), buff5[i].getGroup_name()));
            }
                
            return new ResData(true, gson.toJson(contact)); 
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getContact(FilterContact, int, int). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getContact(FilterContact, int, int). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getContact(FilterContact, int, int). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getContact(FilterContact, int, int). Exception", name_class);
        }
    }

    public ResData getViewContacts (FilterContact item) throws MyNotAuthorizedException {
        try {
            ViewContact[] contacts = gson.fromJson(Select("SELECT * FROM Contacts " + SqlFilter(item, item.getLimit(), item.getOffset()) + ";", connection_with_db, gson), ViewContact[].class);
            for(int j = 0; j < contacts.length;j++) {
                ContactEmail[] buff2 = gson.fromJson(Select("SELECT * FROM ContactEmails WHERE user_id = " + item.getUser_id() + " AND contact_id = " + contacts[j].getId()+ ";", connection_with_db, gson),ContactEmail[].class);
                for(int i = 0; i < buff2.length; i++) {
                    contacts[j].getEmails().add(buff2[i].getEmail());
                }
            }
            
            for(int j = 0; j < contacts.length;j++) {
                ContactPhone[] buff2 = gson.fromJson(Select("SELECT * FROM ContactPhones WHERE user_id = " + item.getUser_id() + " AND contact_id = " + contacts[j].getId()+ ";", connection_with_db, gson),ContactPhone[].class);
                for(int i = 0; i < buff2.length; i++) {
                    contacts[j].getPhones().add(buff2[i].getPhone());
                }
            }
            
            for(int j = 0; j < contacts.length;j++) {
                ContactSocial[] buff2 = gson.fromJson(Select("SELECT * FROM ContactSocial WHERE user_id = '" + item.getUser_id() + "' AND contact_id = " + contacts[j].getId()+ ";", connection_with_db, gson), ContactSocial[].class);
                for(int i = 0; i < buff2.length; i++) {
                    contacts[j].getSocial().add(buff2[i].getSocial());
                }
            }
            for(int i = 0; i < contacts.length; i++) {
                GroupContacts[] buff5 = gson.fromJson(Select("SELECT c.id, c.group_name FROM Contacts a, Contacts_to_GroupContacts b, GroupsContacts c WHERE c.user_id = '" + contacts[i].getUser_id() + "' AND c.id = b.group_contacts_id AND a.id = b.contacts_id  and a.id = " + contacts[i].getId() + ";", connection_with_db, gson), GroupContacts[].class);
                for(int j = 0; j < buff5.length; j++) {
                    contacts[i].getGroups_contacts().add(new ShortGroupContacts(buff5[j].getId(), buff5[j].getGroup_name()));
                }
            }
            return new ResData(true, gson.toJson(contacts));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getViewContacts(FilterContact, int, int). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getViewContacts(FilterContact, int, int). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getViewContacts(FilterContact, int, int). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getViewContacts(FilterContact, int, int). Exception", name_class);
        }    
    }
    
    public ResData getCountContacts(FilterContact item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Select("SELECT COUNT(*) as result FROM Contacts " + SqlFilter(item, 0, 0) + ";", connection_with_db, gson));
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "getCountContacts(FilterContact). Exeption In parse", name_class);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "getCountContacts(FilterContact). Sql all select", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "getCountContacts(FilterContact). Sql select", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "getCountContacts(FilterContact). Exception", name_class);
        }
    }
    
    //Add
    
    public ResData AddContact(ContactModel item) throws MyNotAuthorizedException {
        String result = "", Bffresult = "";
        try {
            ArrayList<Integer> group_add = new ArrayList<Integer>();
            ArrayList<Integer> group_add2 = new ArrayList<Integer>();
            result = Add(item, "Contacts", connection_with_db, gson, true);
            SQLResult res = gson.fromJson(result, SQLResult.class);
            for(int i = 0; i < item.getEmails().size(); i++) {
                Bffresult = Add(new ContactEmail(item.getUser_id(), res.getId(), item.getEmails().get(i)), "ContactEmails", connection_with_db, gson, true);
            }
            
            for(int i = 0; i < item.getPhones().size(); i++) {
                Bffresult = Add(new ContactPhone(item.getUser_id(), res.getId(), item.getPhones().get(i)), "ContactPhones", connection_with_db, gson, true);
            }
                        
            for(int i = 0; i < item.getSocial().size(); i++) {
                Bffresult = Add(new ContactSocial(item.getUser_id(), res.getId(), item.getSocial().get(i)), "ContactSocial", connection_with_db, gson, true);
            }
            
            for(int i = 0; i < item.getGroups_contacts().size();i++) {
              group_add.add(item.getGroups_contacts().get(i).getId());
              TypeGroupContacts buff = gson.fromJson(Select("select default_group_id from TypesGroupContacts where id = (select type_id from GroupsContacts where id = " + item.getGroups_contacts().get(i).getId() + " and user_id = " + item.getUser_id() + ") and user_id = " + item.getUser_id() + ";", connection_with_db, gson), TypeGroupContacts[].class)[0];
              group_add.add(buff.getDefault_group_id());
            }
            
            TypeGroupContacts buff = gson.fromJson(Select("select id from GroupsContacts where type_id = (select type_group_contacts_id from UserSystemInfo where user_id = " + item.getUser_id() + ") and user_id = " + item.getUser_id() + ";", connection_with_db, gson), TypeGroupContacts[].class)[0];
            group_add.add(buff.getId());
            
            for(int i = 0; i < group_add.size();i++) {
                boolean prov = true;
                for(int j = 0; j < group_add2.size();j++) {
                    if(group_add2.get(j).compareTo(group_add.get(i)) == 0) {
                        prov = false;
                        break;
                    }
                }
                if(prov)
                    group_add2.add(group_add.get(i));
            }

            for(int i = 0; i < group_add2.size();i++) {
                Bffresult = Add(new ContactsToGroupContact(item.getUser_id(), res.getId(), group_add2.get(i)), "Contacts_to_GroupContacts ", connection_with_db, gson, true);
            }

            return new ResData(true, result);
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "AddContact(ContactModel). Sql add", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "AddContact(ContactModel). Exeption in parse", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "AddContact(ContactModel). MyTypeException", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "AddContact(ContactModel). Exception", name_class);
        }
    }
    
    //Delete
    
    public ResData DeleteContact(ViewContact item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Delete(item, "Contacts", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "DeleteContact(ViewContact). Sql delete", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "DeleteContact(ViewContact). Exception", name_class);
        } 
    }
    
     //Update
    
    public ResData UpdateContact(ContactModel item) throws MyNotAuthorizedException {
        try {
            ResData res = new ResData();
            res.setRes(Update(item, "Contacts", connection_with_db, gson));
            ArrayList<ShortGroupContacts> db_list = new ArrayList<ShortGroupContacts>();
            ArrayList<ShortGroupContacts> item_list = item.getGroups_contacts();
            GroupContacts[] buff;
            try {
                buff = gson.fromJson(Select("SELECT c.id, c.group_name FROM Contacts a, Contacts_to_GroupContacts b, GroupsContacts c WHERE c.user_id = '" + item.getUser_id() + "' AND c.id = b.group_contacts_id AND a.id = b.contacts_id  and a.id = " + item.getId() + ";", connection_with_db, gson), GroupContacts[].class);
                for(int i = 0; i < buff.length; i++) {
                    db_list.add(new ShortGroupContacts(buff[i].getId(), buff[i].getGroup_name()));
                }
                
                for(int i = 0; i < db_list.size();i++) {
                    for(int j = 0; j < item_list.size();j++) {
                        if(db_list.get(i).getId() == item_list.get(j).getId()) {
                            db_list.remove(i);
                            item_list.remove(j);
                            i = -1;
                            break;
                        }
                    }
                }
                res = Add_contact_in_new_group(item, item_list);
                if(!res.getIs_success()) {
                    return new ResData(false, res.getRes());
                }
                
                res = Delete_contact_in_old_group(item, db_list);
                if(!res.getIs_success()) {
                    return new ResData(false, res.getRes());
                } 
                
                ExecuteQuery("DELETE FROM ContactPhones  where contact_id = " + item.getId() + " and user_id = " + item.getUser_id() + ";", connection_with_db, gson);
                ExecuteQuery("DELETE FROM ContactEmails  where contact_id = " + item.getId() + " and user_id = " + item.getUser_id() + ";", connection_with_db, gson);
                
                for(int i = 0; i < item.getPhones().size();i++) {
                    Add(new ContactPhone(item.getUser_id(), item.getId(), item.getPhones().get(i)), "ContactPhones", connection_with_db, gson, true);
                }
                
                for(int i = 0; i < item.getEmails().size();i++) {
                    Add(new ContactEmail(item.getUser_id(), item.getId(), item.getEmails().get(i)), "ContactEmails", connection_with_db, gson, true);
                }
            } 
            catch (MyTypeException ex) {
                return log_manager.CallError(ex.toString(), "UpdateContact(ContactModel). MyTypeException", name_class);
            } 
            return new ResData(true, res.getRes());
            
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateContact(ContactModel). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateContact(ContactModel). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateContact(ContactModel). Exception", name_class);
        }
    }
    
    public ResData UpdateStatusContact(ContactModel item) throws MyNotAuthorizedException {
        try {
            return new ResData(true, Update(item, "Contacts", connection_with_db, gson));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "UpdateContact(ContactModel). Sql update", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "UpdateContact(ContactModel). Exception In parse", name_class);
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "UpdateContact(ContactModel). Exception", name_class);
        }
    }
    
    private ResData Add_contact_in_new_group(ContactModel item, ArrayList<ShortGroupContacts> add_list) {
        int size = add_list.size();
        for(int i = 0; i < size; i ++) {
            try {
                TypeGroupContacts buff = gson.fromJson(Select("select default_group_id from TypesGroupContacts where id = (select type_id from GroupsContacts where id = " + add_list.get(i).getId() + ") and user_id = " + item.getUser_id() + ";", connection_with_db, gson), TypeGroupContacts[].class)[0];
                boolean prov = true;
                for(int j = 0; j < add_list.size();j++) {
                    if(add_list.get(j).getId() == buff.getDefault_group_id()) {
                        prov = false;
                        break;
                    }
                }
                if(prov) {
                    add_list.add(new ShortGroupContacts(buff.getDefault_group_id(), ""));
                }
            } catch (SQLException ex) {
                return log_manager.CallError(ex.toString(), "Add_contact_in_new_group(ContactModel, ArrayList<ShortGroupContacts>). SQLException", name_class);
            } catch (MyTypeException ex) {
                return log_manager.CallError(ex.toString(), "Add_contact_in_new_group(ContactModel, ArrayList<ShortGroupContacts>). MyTypeException", name_class);
            } catch (Exception ex) {
                return log_manager.CallError(ex.toString(), "Add_contact_in_new_group(ContactModel, ArrayList<ShortGroupContacts>). Exception", name_class);
            }
        }
        for(int i = 0; i < add_list.size();i++) {
            try {
                Add(new ContactsToGroupContact(item.getUser_id(),item.getId(),add_list.get(i).getId()),"Contacts_to_GroupContacts ",connection_with_db, gson, true);
            } catch (SQLException ex) {
                if(ex.getErrorCode() != 1062) {
                    return log_manager.CallError(ex.toString(), "Add_contact_in_new_group(ContactModel, ArrayList<ShortGroupContacts>). SQLException 2 " + ex.getErrorCode(), name_class);    
                }
            } catch (MyParseException ex) {
                return log_manager.CallError(ex.toString(), "Add_contact_in_new_group(ContactModel, ArrayList<ShortGroupContacts>). MyParseException 2", name_class);
            } catch (MyNotAuthorizedException ex) {
                return log_manager.CallError(ex.toString(), "Add_contact_in_new_group(ContactModel, ArrayList<ShortGroupContacts>). MyNotAuthorizedException ", name_class);
            } catch (Exception ex) {
                return log_manager.CallError(ex.toString(), "Add_contact_in_new_group(ContactModel, ArrayList<ShortGroupContacts>). Exception 2", name_class);
            }
        }
        
        return new ResData(true, gson.toJson(new SQLResult("Ok"), SQLResult.class));
    }
    
    private ResData Delete_contact_in_old_group(ContactModel item, ArrayList<ShortGroupContacts> delete_list) {
        try {
            int main_group_def = -1;
            
            UserSystemInfo buff = gson.fromJson(Select("select id from GroupsContacts where type_id = (select type_group_contacts_id from UserSystemInfo where user_id = " + item.getUser_id() + ");", connection_with_db, gson), UserSystemInfo[].class)[0];
            main_group_def = buff.getId();
                
            for(int i = 0; i < delete_list.size();i++) {
                if(delete_list.get(i).getId() == main_group_def) {
                    ContactsToGroupContact buff3[] = gson.fromJson(Select("select * from Contacts_to_GroupContacts where contacts_id = " + item.getId() + ";", connection_with_db, gson), ContactsToGroupContact[].class);
                    if(buff3.length == 1)
                        ExecuteQuery("delete from Contacts_to_GroupContacts where contacts_id = " + item.getId() + " and group_contacts_id = " + delete_list.get(i).getId() + ";", connection_with_db, gson);
                }
                else {
                    TypeGroupContacts buff2[] = gson.fromJson(Select("select id, default_group_id from TypesGroupContacts where user_id = " + item.getUser_id() + " and default_group_id = " + delete_list.get(i).getId() + ";", connection_with_db, gson), TypeGroupContacts[].class);
                    if(buff2.length == 1) {
                        GroupContacts buff3[] = gson.fromJson(Select("select * from Contacts_to_GroupContacts where contacts_id = " + item.getId()
                                + " and group_contacts_id in (select id from GroupsContacts where type_id  = " + buff2[0].getId() + ");", connection_with_db, gson), GroupContacts[].class);
                        if(buff3.length == 1) 
                            ExecuteQuery("delete from Contacts_to_GroupContacts where contacts_id = " + item.getId() + " and group_contacts_id = " + delete_list.get(i).getId() + ";", connection_with_db, gson);
                    }
                    else {
                        ExecuteQuery("delete from Contacts_to_GroupContacts where contacts_id = " + item.getId() + " and group_contacts_id = " + delete_list.get(i).getId() + ";", connection_with_db, gson);
                    }
                }
            }
            return new ResData(true, gson.toJson(new SQLResult("Ok"), SQLResult.class));
        } catch (SQLException ex) {
            return log_manager.CallError(ex.toString(), "Delete_contact_in_old_group(ContactModel, ArrayList<ShortGroupContacts>). SQLException", name_class);
        } catch (MyTypeException ex) {
            return log_manager.CallError(ex.toString(), "Delete_contact_in_old_group(ContactModel, ArrayList<ShortGroupContacts>). MyTypeException", name_class);
        } catch (MyParseException ex) {
            return log_manager.CallError(ex.toString(), "Delete_contact_in_old_group(ContactModel, ArrayList<ShortGroupContacts>). MyParseException", name_class);
        } catch (MyNotAuthorizedException ex) {
            return log_manager.CallError(ex.toString(), "Delete_contact_in_old_group(ContactModel, ArrayList<ShortGroupContacts>). MyNotAuthorizedException", name_class);
        } catch (Exception ex){
            return log_manager.CallError(ex.toString(), "Delete_contact_in_old_group(ContactModel, ArrayList<ShortGroupContacts>). Exception", name_class);
        }
    }
    
    public ResData Contact_select_with_count(FilterContact item){
        try{
            ResData res, answer;
            
            answer = getViewContacts(item);
            if(!answer.getIs_success())
                return answer;
            
            res = getCountContacts(item);
            if(!res.getIs_success())
                return res;
            
            answer.setRes("[" + answer.getRes() + "," + res.getRes().substring(1, res.getRes().length() - 1) + "]");
            answer.setIs_success(res.getIs_success());
            
            return answer;
        }
        catch(MyNotAuthorizedException ex){
            return log_manager.CallError(ex.toString(), "Contact_select_with_count(FilterContact item). MyNotAuthorizedException", name_class);
        } catch (Exception ex){
            return log_manager.CallError(ex.toString(), "Contact_select_with_count(FilterContact item). Exception", name_class);
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
