package SQLClass.VkError;

import com.google.gson.JsonObject;


public class VkError {
    String err;
    
    public VkError(String err) {
        this.err = err;
    }
    
    public String getErr() {
        return err;
    }

    public void setErr(String err) {
        this.err = err;
    }
    
    
    
}
