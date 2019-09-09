/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SQLClass.VkError;

public class VkErrorAdvanced {
    int code;
    String error_msg;
    String request_params = "";
    
    public VkErrorAdvanced(int code, String error_msg, String request_params){
        this.code = code;
        this.error_msg = error_msg;
        this.request_params = request_params;
    }

    public VkErrorAdvanced(int code, String error_msg){
        this.code = code;
        this.error_msg = error_msg;
    }
    
    public int getCode() {
        return code;
    }

    public void setCode(int code) {
        this.code = code;
    }

    public String getError_msg() {
        return error_msg;
    }

    public void setError_msg(String error_msg) {
        this.error_msg = error_msg;
    }

    public String getRequest_params() {
        return request_params;
    }

    public void setRequest_params(String request_params) {
        this.request_params = request_params;
    }
    
    
}
