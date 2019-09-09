package SQLClass.Contact;

import MyInterface.ModelInterface;
import java.sql.Timestamp;
import java.util.ArrayList;

public class ViewContact implements ModelInterface {
    int id = -1; 
    int user_id = -1;
    String surname = "";
    String name = "";
    String patronymic = "";
    String contact_status = "";
    Timestamp date_create;
    
    ArrayList<ShortGroupContacts> groups_contacts;
    ArrayList<String> social;
    ArrayList<String> emails;
    ArrayList<String> phones;
    
    public ViewContact(int id, int user_id, String surname, String name, String patronymic, String contact_status, 
            Timestamp date_create, ArrayList<ShortGroupContacts> groups_contacts, ArrayList<String> social, ArrayList<String> emails, ArrayList<String> phones) {
        this.id = id; 
        this.user_id = user_id;
        this.surname = surname;
        this.name = name;
        this.patronymic = patronymic;
        this.contact_status = contact_status;
        this.date_create = date_create;
        this.groups_contacts = groups_contacts;
        this.social = social;
        this.emails = emails;
        this.phones = phones;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), surname, name, patronymic, contact_status, date_create.toString()};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id),
            String.valueOf(user_id), 
            "'" + surname.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + patronymic.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + contact_status.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + date_create.toString() + "'",
        };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "surname", "name", "patronymic", "contact_status", "date_create"}; 
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public String getSurname() {
        return surname;
    }

    public void setSurname(String surname) {
        this.surname = surname;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getPatronymic() {
        return patronymic;
    }

    public void setPatronymic(String patronymic) {
        this.patronymic = patronymic;
    }

    public String getContact_status() {
        return contact_status;
    }

    public void setContact_status(String contact_status) {
        this.contact_status = contact_status;
    }

    public Timestamp getDate_create() {
        return date_create;
    }

    public void setDate_create(Timestamp date_create) {
        this.date_create = date_create;
    }

    public ArrayList<ShortGroupContacts> getGroups_contacts() {
        return groups_contacts;
    }

    public void setGroups_contacts(ArrayList<ShortGroupContacts> groups_contacts) {
        this.groups_contacts = groups_contacts;
    }

    public ArrayList<String> getSocial() {
        return social;
    }

    public void setSocial(ArrayList<String> social) {
        this.social = social;
    }

    public ArrayList<String> getEmails() {
        return emails;
    }

    public void setEmails(ArrayList<String> emails) {
        this.emails = emails;
    }

    public ArrayList<String> getPhones() {
        return phones;
    }

    public void setPhones(ArrayList<String> phones) {
        this.phones = phones;
    }
}
