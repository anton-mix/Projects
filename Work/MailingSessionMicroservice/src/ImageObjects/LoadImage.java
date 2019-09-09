package ImageObjects;

import java.util.ArrayList;

public class LoadImage {
    String user_id = "";
    int microservice_id = -1;
    int image_id = -1;
    ArrayList<String> part_path;
   
    public LoadImage(String user_id, int microservice_id, int image_id, ArrayList<String> part_path) {
        this.user_id = user_id;
        this.microservice_id = microservice_id;
        this.image_id = image_id;
        this.part_path = part_path;
    }
    
    public LoadImage(String user_id, int microservice_id, int image_id) {
        this.user_id = user_id;
        this.microservice_id = microservice_id;
        this.image_id = image_id;
        part_path = new ArrayList<String>();
    }

    public String getUser_id() {
        return user_id;
    }

    public void setUser_id(String user_id) {
        this.user_id = user_id;
    }

    public int getMicroservice_id() {
        return microservice_id;
    }

    public void setMicroservice_id(int microservice_id) {
        this.microservice_id = microservice_id;
    }

    public int getImage_id() {
        return image_id;
    }

    public void setImage_id(int image_id) {
        this.image_id = image_id;
    }

    public ArrayList<String> getPart_path() {
        return part_path;
    }

    public void setPart_path(ArrayList<String> part_path) {
        this.part_path = part_path;
    }
}
