package SQLClass.GroupObject;


public class GroupObject {
    int id;
    String name;
    String screen_name;
    int is_closed;
    String type;
    int is_admin;
    int admin_level;
    int is_member;
    int is_advertiser;
    String photo_50;
    String photo_100;
    String photo_200;
    
    
    public GroupObject(){
         
    }
    
    public GroupObject(int id, String name, String screen_name, int is_closed, String type,
        int is_admin, int admin_level, int is_member, int is_advertiser, String photo_50, String photo_100, String photo_200){
            this.id = id;
            this.name = name;
            this.screen_name = screen_name;
            this.is_closed = is_closed;
            this.type = type;
            this.is_admin = is_admin;
            this.admin_level = admin_level;
            this.is_member = is_member;
            this.is_advertiser = is_advertiser;
            this.photo_50 = photo_50;
            this.photo_100 = photo_100;
            this.photo_200 = photo_200;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getScreen_name() {
        return screen_name;
    }

    public void setScreen_name(String screen_name) {
        this.screen_name = screen_name;
    }

    public int getIs_closed() {
        return is_closed;
    }

    public void setIs_closed(int is_closed) {
        this.is_closed = is_closed;
    }

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }

    public int getIs_admin() {
        return is_admin;
    }

    public void setIs_admin(int is_admin) {
        this.is_admin = is_admin;
    }

    public int getAdmin_level() {
        return admin_level;
    }

    public void setAdmin_level(int admin_level) {
        this.admin_level = admin_level;
    }

    public int getIs_member() {
        return is_member;
    }

    public void setIs_member(int is_member) {
        this.is_member = is_member;
    }

    public int getIs_advertiser() {
        return is_advertiser;
    }

    public void setIs_advertiser(int is_advertiser) {
        this.is_advertiser = is_advertiser;
    }

    public String getPhoto_50() {
        return photo_50;
    }

    public void setPhoto_50(String photo_50) {
        this.photo_50 = photo_50;
    }

    public String getPhoto_100() {
        return photo_100;
    }

    public void setPhoto_100(String photo_100) {
        this.photo_100 = photo_100;
    }

    public String getPhoto_200() {
        return photo_200;
    }

    public void setPhoto_200(String photo_200) {
        this.photo_200 = photo_200;
    }
    
    
}
