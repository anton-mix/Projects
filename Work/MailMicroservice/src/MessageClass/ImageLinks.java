package MessageClass;

import MyInterface.ModelInterface;


public class ImageLinks implements ModelInterface{
    int id = -1;
    int user_id = -1;
    int message_id = -1;
    String image_address = "";
    String image_name = "";
    
    public ImageLinks(int id, int user_id, int message_id, String image_address, String image_name) {
        this.id = id;
        this.user_id = user_id;
        this.message_id = message_id;
        this.image_address = image_address;
        this.image_name = image_name;
    }
    
    public ImageLinks(int user_id, int message_id, String image_address, String image_name) {
        this.user_id = user_id;
        this.message_id = message_id;
        this.image_address = image_address;
        this.image_name = image_name;
    }
    
    public ImageLinks(String image_address, String image_name) {
        this.image_address = image_address;
        this.image_name = image_name;
    }
        
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(message_id), 
            image_address, image_name};
    }
    
    @Override
     public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id),
            String.valueOf(message_id), 
            "'" + image_address.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + image_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"
        };
    }
     
     @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "message_id", "image_address", "image_name"}; 
    }

    @Override
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

    public int getMessage_id() {
        return message_id;
    }

    public void setMessage_id(int message_id) {
        this.message_id = message_id;
    }

    public String getImage_address() {
        return image_address;
    }

    public void setImage_address(String image_address) {
        this.image_address = image_address;
    }

    public String getImage_name() {
        return image_name;
    }

    public void setImage_name(String image_name) {
        this.image_name = image_name;
    }
     
    
}