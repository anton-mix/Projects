package SQLClass.Contact;

public class ShortContact {
    int id = -1;
    String contact_name = "";
    
    public ShortContact() {
        id = -1;
        contact_name = "";
    }
    
    public ShortContact(int id, String contact_name) {
        this.id = id;
        this.contact_name = contact_name;
    }

    public String getContact_name() {
        return contact_name;
    }

    public void setContact_name(String contact_name) {
        this.contact_name = contact_name;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }
}
