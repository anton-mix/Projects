package SQLClass.Contact;

import MyInterface.ModelInterface;
import java.sql.Timestamp;
import java.util.ArrayList;

public class FilterContact implements ModelInterface{
    int user_id = -1;
    String city = "";
    String client_info = "";
    String tags = "";
    Timestamp date_create;
    String query = "";
    ArrayList<ShortGroupContacts> groups_contacts;
    ArrayList<String> social;
    ArrayList<String> emails;
    ArrayList<String> phones;
    int limit = 0;
    int offset = 0;  
    
    public FilterContact(int user_id, String city, String client_info, String tags, Timestamp date_create, String query, 
            ArrayList<ShortGroupContacts> groups_contacts ,ArrayList<String> social, ArrayList<String> emails, ArrayList<String> phones, int limit, int offset) {
        this.user_id = user_id;
        this.city = city;
        this.client_info = client_info;
        this.tags = tags;
        this.date_create = date_create;
        this.query = query;
        this.groups_contacts = groups_contacts;
        this.social = social;
        this.emails = emails;
        this.phones = phones;
        this.limit = limit;
        this.offset = offset;
    }
    
    public String getFullText() {
        return "name, surname, patronymic, contact_status";
    }
    
    public String[] getArray() {
        return new String[] {String.valueOf(user_id), city, client_info, tags, date_create.toString(), query, String.valueOf(limit), String.valueOf(offset)};
    }
    
    public String[] getSqlArray() {
        return new String[] {
            String.valueOf(user_id), 
            "'" + city.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + client_info.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + tags.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + date_create.toString() + "'"
            };
    }
     
    @Override
    public String[] getKeyArray() {
        return new String[] {"user_id", "city", "client_info", "tags", "date_create", "query", "limit", "offset"}; 
    }

    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public String getCity() {
        return city;
    }

    public void setCity(String city) {
        this.city = city;
    }

    public String getClient_info() {
        return client_info;
    }

    public void setClient_info(String client_info) {
        this.client_info = client_info;
    }

    public String getTags() {
        return tags;
    }

    public void setTags(String tags) {
        this.tags = tags;
    }

    public Timestamp getDate_create() {
        return date_create;
    }

    public void setDate_create(Timestamp date_create) {
        this.date_create = date_create;
    }

    @Override
    public String getQuery() {
        return query;
    }

    public void setQuery(String query) {
        this.query = query;
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

    @Override
    public int getLimit() {
        return limit;
    }

    public void setLimit(int limit) {
        this.limit = limit;
    }

    @Override
    public int getOffset() {
        return offset;
    }

    public void setOffset(int offset) {
        this.offset = offset;
    }
    
    
}
