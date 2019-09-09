package SQLClass;

public class ResData {
    boolean is_success = false;
    String res;
    
    public ResData(){
        this.is_success = false;
        this.res = "";
    }

    public ResData(boolean is_success,  String res){
        this.is_success = is_success;
        this.res = res;
    }
    
    public boolean getIs_success() {
        return is_success;
    }

    public void setIs_success(boolean is_success) {
        this.is_success = is_success;
    }

    public String getRes() {
        return res;
    }

    public void setRes(String res) {
        this.res = res;
    }
    
    
}
