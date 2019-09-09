package SQLClass.Contact;

public class ViewTypeGroupContacts {
    int id = -1;
    int user_id = -1; 
    String name = ""; 
    
    public ViewTypeGroupContacts (int id, int user_id, String name)  {
        this.id = id;
        this.user_id = user_id;
        this.name = name;   
    }
    
    public ViewTypeGroupContacts (int user_id, String name)  {
        this.user_id = user_id;
        this.name = name;   
    }
    
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), name};
    }
    
    public String[] getSqlArray() {
            return new String[] {String.valueOf(id),
                String.valueOf(user_id),
                "'" + name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
        };
    }
    
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "name"}; 
    }
    
    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
}
