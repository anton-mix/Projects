package MyThreads;

import Managers.MessageManager;
import MyEnum.RestRamEnum;
import RestObject.RestObject;
import SQLClass.SQLResult;
import com.google.gson.Gson;
import com.google.gson.JsonObject;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;

public class WriteThread implements Runnable {
    private String ErrorString = "";
    private String EventString = "";
    long time_refresh = 86400000;
    long cur_time = time_refresh;
    private ArrayList<String> ErrorAll;
    private ArrayList<String> EventAll;
    private ArrayList<String> Hash_map_session;
    private ArrayList<String> User_auth_session;
    File file;
    static DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd");
    MessageManager message_manager;
    Gson gson;
    SQLResult res;
    
    public WriteThread(MessageManager message_manager, Gson gson) {
        dateFormat = new SimpleDateFormat("yyyy.MM.dd");
        System.out.println("Start MyThread");
        this.message_manager = message_manager;
        this.gson = gson;
        Hash_map_session = new ArrayList<String>();
        User_auth_session = new ArrayList<String>();
        EventAll = new ArrayList<String>();
        ErrorAll = new ArrayList<String>();
    }
    
    @Override
    public void run() {
        while (true) {
            try {
                Date date = new Date();
                
                file = new File("logs/" + dateFormat.format(date));
                if(!file.exists()) {
                    file.mkdir();
                }
                System.out.println("Writing....");
                if(EventAll.size() != 0) {
                    doWork(EventAll, 1);
                    //EventAll = null;
                }
                if(ErrorAll.size() != 0) {
                    doWork(ErrorAll, 2);
                    //ErrorAll = null;
                }
                
                try {
                    res = gson.fromJson(message_manager.call(gson.toJson(new RestObject(new JsonObject(), "", RestRamEnum.Get_session_keys.ordinal())), "ram_queue"), SQLResult.class);
                
                } catch (IOException ex) {
                    System.out.println("Error " + ex.toString());
                    Logger.getLogger(WriteThread.class.getName()).log(Level.SEVERE, null, ex);
                }
                
                if(!"".equals(res.getResult())) {
                    Hash_map_session.add(res.getResult());
                    doWork(Hash_map_session, 3);
                }
                //Hash_map_session.remove(0);
                res.setResult("");
                
                try {
                    res = gson.fromJson(message_manager.call(gson.toJson(new RestObject(new JsonObject(), "", RestRamEnum.Get_auth_message_keys.ordinal())), "ram_queue"), SQLResult.class);
                } catch (IOException ex) {
                    System.out.println("Error " + ex.toString());
                    Logger.getLogger(WriteThread.class.getName()).log(Level.SEVERE, null, ex);
                }
                
                if(!"".equals(res.getResult())) {
                    User_auth_session.add(res.getResult());
                    doWork(User_auth_session, 4);
                }
                
                //password messages
                res.setResult("");
                
                try {
                    res = gson.fromJson(message_manager.call(gson.toJson(new RestObject(new JsonObject(), "", RestRamEnum.Get_change_pass_message_keys.ordinal())), "ram_queue"), SQLResult.class);
                } catch (IOException ex) {
                    System.out.println("Error " + ex.toString());
                    Logger.getLogger(WriteThread.class.getName()).log(Level.SEVERE, null, ex);
                }
                
                if(!"".equals(res.getResult())) {
                    User_auth_session.add(res.getResult());
                    doWork(User_auth_session, 5);
                }
                
                
                //User_auth_session.remove(0);
                res.setResult("");
                
                if(cur_time >= time_refresh) {
                    try {
                        res = gson.fromJson(message_manager.call(gson.toJson(new RestObject(new JsonObject(), "", RestRamEnum.Refres_queues.ordinal())), "ram_queue"), SQLResult.class);
                    } catch (IOException ex) {
                        System.out.println("Error " + ex.toString());
                        Logger.getLogger(WriteThread.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    cur_time = 0;
                }
                
                Thread.sleep(30000);
                cur_time = cur_time + 30000;
            } catch (InterruptedException ex) {
                System.out.println("Error " + ex.toString());
                System.out.println("Thread sleep " + ex.toString());
            }
        }
    }
    
    private static void doWork(ArrayList<String> message, int type) {
        FileWriter writer = null;
        File file;
        String buff = "";
        Date date = new Date();
        while (message.size() > 0){
            buff = message.get(0) + "\n";
            for(int i = 1; i < message.size(); i++){
                if((message.get(0).split(" ")[10].split(",")[0]).equals(message.get(i).split(" ")[10].split(",")[0])){
                    buff = buff + message.get(i) + "\n";
                    message.remove(i);
                    i--;
                }
            }
            
            if (type <= 2){
                file = new File("logs/" + dateFormat.format(date) + "/" + (message.get(0).split(" ")[10].split(",")[0]));   
            }
            else{
                file = new File("logs/" + dateFormat.format(date) + "/UserSessionMicroservice");
            }
            
            if(!file.exists()) {
                file.mkdir();
            }
      
            message.remove(0);
                        
            try {
                switch (type) {
                    case 1:
                        writer = new FileWriter(file.getPath() + "/event.txt",true);
                        break;
                    case 2:
                        writer = new FileWriter(file.getPath() + "/error.txt", true);
                        break;
                    case 3:
                        writer = new FileWriter(file.getPath() + "/session.txt", false);
                        break;
                    case 4:
                        writer = new FileWriter(file.getPath() + "/user_auth.txt", false);
                        break;
                    case 5:
                        writer = new FileWriter(file.getPath() + "/password_messages.txt", false);
                        break;
                    default:
                        System.out.println("Error type in writer");
                        break;
                }
                writer.write(buff);
                writer.flush();
                writer.close();
                buff = "";
            } catch (IOException ex) {
                System.out.println("error: " + ex.getMessage());
            }
        }
    }
    
    public String getErrorString() {
        return ErrorString;
    }

    public void setErrorString(String ErrorString) {
        this.ErrorString = ErrorString;
        ErrorAll.add(ErrorString);
    }

    public String getEventString() {
        return EventString;
    }
    
    public void setEventString(String EventString) {
        this.EventString = EventString;
        EventAll.add(EventString);
    }
    
    
}
