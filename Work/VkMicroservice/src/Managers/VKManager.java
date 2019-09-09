package Managers;

import SQLClass.AnswerGroupListObject.AnswerGroupListObject;
import SQLClass.AnswerGroupListObject.AnswerGroupListObjectConverter;
import SQLClass.AnswerMemberGroupObject.AnswerMemberGroupObject;
import SQLClass.AnswerMemberGroupObject.AnswerMemberGroupObjectConverter;
import SQLClass.VkError.VkErrorConverter;
import SQLClass.GroupObject.GroupObject;
import SQLClass.GroupObject.GroupObjectConverter;
import SQLClass.GlobalVariables;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.VkError.VkError;
import SQLClass.VkError.VkErrorAdvanced;
import SQLClass.VkError.VkErrorAdvancedConverter;
import SQLClass.VKAuthClass.VKAuthClass;
import SQLClass.VKSendingClass.VKSendingClass;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.util.ArrayList;
import java.util.Random;
import TestSender.TestSender;

public class VKManager  {

    public static final String REDIRECT_URL = "https://oauth.vk.com/blank.html";
    public static final String VK_AUTH_URL = "https://oauth.vk.com/authorize?client_id=6793765&display=page&redirect_uri=https://oauth.vk.com/blank.html&scope=266242&response_type=token&v=5.92&state=123456"; //TODO!!!
    public static final String VK_AUTH_URL_GROUP = "https://oauth.vk.com/authorize?client_id=6793765&group_ids=171469935&display=page&redirect_uri=https://oauth.vk.com/blank.html&scope=4096&response_type=token&v=5.92";
    public static String tokenUrl;
    public static String answer;
    public static String user_id;
    public static String tokenURL;
    public static String tokenGrURL;
    GlobalVariables global_var;
    LogManager log_manager;
    Gson gson;
    static int  min = 1000000000;  
    String class_name = "VKManager";
    String ErrorHttp = "";
    
    public VKManager(GlobalVariables global_var, LogManager log_manager) {
        this.log_manager = log_manager;
        this.global_var = global_var;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(AnswerGroupListObject.class, new AnswerGroupListObjectConverter());  
        builder.registerTypeAdapter(GroupObject.class, new GroupObjectConverter()); 
        builder.registerTypeAdapter(AnswerMemberGroupObject.class, new AnswerMemberGroupObjectConverter()); 
        builder.registerTypeAdapter(VkError.class, new VkErrorConverter()); 
        builder.registerTypeAdapter(VkErrorAdvanced.class, new VkErrorAdvancedConverter()); 
        gson = builder.create(); 
    }
    
    public ResData getUserGroups(VKAuthClass vk_auth_class){
        GroupObject[] groups;
        AnswerGroupListObject answ_obj;
        String groups_list;
        VkError err;
        try {  
            groups_list = TestSender.sendPost("https://api.vk.com/method/groups.get?user_id=" + vk_auth_class.getId() + "&extended=1&filter=moder&v=5.92&access_token=" + vk_auth_class.getToken(), null);
            err = gson.fromJson(groups_list, VkError.class);
            if (!err.getErr().equals(""))
                return new ResData(false, gson.toJson(new SQLResult(err.getErr())));
            else{
                groups_list = groups_list.split(":", 2)[1];
                groups_list = groups_list.substring(0, groups_list.length() - 1);          
                answ_obj = gson.fromJson(groups_list, AnswerGroupListObject.class);
                groups = gson.fromJson(answ_obj.getJson_message(), GroupObject[].class);
                log_manager.CallEvent("getUserGroups Success", class_name);
                String res = gson.toJson(groups, GroupObject[].class);
                return new ResData(true, res);
            }
        }
        catch (Exception ex) {
            return log_manager.CallError(ex.toString(), "Error in getUserGroups", class_name);
        } 
    }
    /*
    public ResData sendMessagesByUser(VKSendingClass vk_sending_class){
        //User ids, text, user token
        Random random = new Random();     
              
        try{
            for(int i = 0; i < vk_sending_class.getUser_ids().size(); i++){
                int randomInt = Math.abs(random.nextInt(1000000000)) + min;
                answer = TestSender.sendPost("https://api.vk.com/method/messages.send?user_id=" + vk_sending_class.getUser_ids().get(i) + "&random_id=" + randomInt + "&message=" + vk_sending_class.getMessage_text() + "&access_token=" + vk_sending_class.getToken() + "&v=5.92", null);
            }
            log_manager.CallEvent("sendMessagesByUser Success", class_name);
        }
        catch(Exception ex){
            return log_manager.CallError(ex.toString(), "Error in sendMessagesByUser", class_name);
        }
        return answer;
    }
    **/
    public ResData sendMessagesByGroup(VKSendingClass vk_sending_class){
        //User ids, text, group token
        vk_sending_class.setGroup_id(vk_sending_class.getGroup_token().split(" ",2)[1]); 
        vk_sending_class.setGroup_token(vk_sending_class.getGroup_token().split(" ",2)[0]);
        Random random = new Random();  
        VkError err;
        AnswerMemberGroupObject members;
        AnswerMemberGroupObject members_in_loop = null;
        try{

            String members_list = TestSender.sendPost("https://api.vk.com/method/groups.getMembers?group_id=" + vk_sending_class.getGroup_id() + "&sort=id_asc&offset=0&v=5.92&access_token=" + vk_sending_class.getGroup_token(), null);
            err = gson.fromJson(members_list, VkError.class);
            System.out.println(err.getErr());
            if (!err.getErr().equals("{}"))
                return new ResData(false, gson.toJson(new SQLResult(err.getErr())));
            else{
                System.out.println(0);
                members_list = members_list.split(":", 2)[1];
                members_list = members_list.substring(0, members_list.length() - 1);
                members = gson.fromJson(members_list, AnswerMemberGroupObject.class);
                for (int i = 0; i < members.getCount(); i += 1000){
                    System.out.println(1);
                    members_list = TestSender.sendPost("https://api.vk.com/method/groups.getMembers?group_id=" + vk_sending_class.getGroup_id() + "&sort=id_asc&offset=" + i + "&v=5.92&access_token=" + vk_sending_class.getGroup_token(), null);
                    err = gson.fromJson(members_list, VkError.class);
                    if (!err.getErr().equals("{}"))
                        return new ResData(false, gson.toJson(new SQLResult(err.getErr())));
                    else{
                        System.out.println(2);
                        members_list = members_list.split(":", 2)[1];
                        members_list = members_list.substring(0, members_list.length() - 1);
                        members_in_loop = gson.fromJson(members_list, AnswerMemberGroupObject.class);
                    }
                }

                for(int i = 0; i < members_in_loop.getIds().size(); i++){
                    int randomInt = Math.abs(random.nextInt(1000000000)) + min;
                    if  (!members_in_loop.getIds().get(i).equals("159212321")){
                        answer = TestSender.sendPost("https://api.vk.com/method/messages.send?user_id=" + members_in_loop.getIds().get(i) + "&random_id=" + randomInt + "&message=" + vk_sending_class.getMessage_text() + "&access_token=" + vk_sending_class.getGroup_token() + "&v=5.92", null);
                    }
                    err = gson.fromJson(answer, VkError.class);
                    if (!err.getErr().equals("{}")){
                        if (gson.fromJson(err.getErr(), VkErrorAdvanced.class).getCode() != 901)
                            return new ResData(false, gson.toJson(new SQLResult(err.getErr())));
                    }
                }
                log_manager.CallEvent("sendMessagesByGroup Success", class_name);

            }
        }
        catch(Exception ex){
            return log_manager.CallError(ex.toString(), "Error in sendMessagesByGroup", class_name);
        }
        return new ResData(true, gson.toJson(new SQLResult("Ok")));
    }
 /**
    public ResData getGroupMembers(VKAuthClass vk_auth_class){
        //group id, user/group token
        AnswerMemberGroupObject members;
        AnswerMemberGroupObject members_in_loop;
        ArrayList<Integer> result = new ArrayList<>();
        try {   
            String members_list = TestSender.sendPost("https://api.vk.com/method/groups.getMembers?group_id=" + vk_auth_class.getId() + "&sort=id_asc&offset=0&v=5.92&access_token=" + vk_auth_class.getToken(), null);
            members_list = members_list.split(":", 2)[1];
            members_list = members_list.substring(0, members_list.length() - 1);
            members = gson.fromJson(members_list, AnswerMemberGroupObject.class);
            for (int i = 0; i < members.getCount(); i += 1000){
                members_list = TestSender.sendPost("https://api.vk.com/method/groups.getMembers?group_id=" + vk_auth_class.getId() + "&sort=id_asc&offset=" + i + "&v=5.92&access_token=" + vk_auth_class.getToken(), null);
                members_list = members_list.split(":", 2)[1];
                members_list = members_list.substring(0, members_list.length() - 1);
                members_in_loop = gson.fromJson(members_list, AnswerMemberGroupObject.class);
                result.addAll(members_in_loop.getIds());
            }
            for(int i = 0; i < vk_sending_class.getUser_ids().size(); i++){
                int randomInt = Math.abs(random.nextInt(1000000000)) + min;
                answer = TestSender.sendPost("https://api.vk.com/method/messages.send?user_id=" + vk_sending_class.getUser_ids().get(i) + "&random_id=" + randomInt + "&message=" + vk_sending_class.getMessage_text() + "&access_token=" + vk_sending_class.getToken() + "&v=5.92", null);
            }
            
            log_manager.CallEvent("getGroupMembers Success", class_name);
            String res = gson.toJson(result, ArrayList[].class);
            return new ResData(true, res);
        }
        catch(Exception ex){
            return log_manager.CallError(ex.toString(), "Error in getGroupMembers", class_name);
        }
    }
    **/
    public String joinUserToGroup(VKAuthClass vk_auth_class){
        //User token + Group Id
        String responce = TestSender.sendPost("https://api.vk.com/method/groups.join?group_id=" + vk_auth_class.getId() + "&access_token="+ vk_auth_class.getToken() +"&v=5.92", null);
        return responce;
    }
    
    public String allowMessagesFromGroup(VKAuthClass vk_auth_class){
        String responce = TestSender.sendPost("https://api.vk.com/method/messages.allowMessagesFromGroup?group_id=" + vk_auth_class.getId() + "&v=5.92&access_token=" + vk_auth_class.getToken(), null);
        return responce;
    }
    
    public String getErrorHttp() {
        return ErrorHttp;
    }

    public void setErrorHttp(String ErrorHttp) {
        this.ErrorHttp = ErrorHttp;
    }
/*
    public static String getToken() {
        String token;
        launch(VKManager.class);
        token = tokenUrl.split("=")[1].split("&")[0];
        return token;
    }

    @Override
    public void start(Stage primaryStage) throws Exception {
        final WebView view = new WebView();
        final WebEngine engine = view.getEngine();
        engine.load(VK_AUTH_URL);
        primaryStage.setScene(new Scene(view));
        primaryStage.show();

        engine.locationProperty().addListener(new ChangeListener<String>() {
            @Override
            public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
                if (newValue.startsWith(REDIRECT_URL)) {
                    tokenUrl = newValue;
                    primaryStage.close();
                }
            }

        });

    }
**/
}
