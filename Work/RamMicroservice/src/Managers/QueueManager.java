package Managers;

import java.io.UnsupportedEncodingException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.sql.Timestamp;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

public class QueueManager {
    Map<String, Object> queue;
    LogManager log_manager;
    
    String name_class = "";
    
    public QueueManager(LogManager log_manager, String name_class) {
        queue = new HashMap<String, Object>();
        this.name_class = name_class;
        this.log_manager = log_manager;
    }

    public Map<String, Object> getQueue() {
        return queue;
    }

    public void setQueue(Map<String, Object> queue) {
        this.queue = queue;
    }
    
    public void AddItem(String key, Object item) {
        queue.put(key, item);
    }
    
    public void DeleteItemFromKey(String key) {
        queue.remove(key);
    }
    
    public void DeleteItemFromItem(Object item) {
        for(int i = 0; i < queue.values().toArray().length;i++) {
            if(queue.values().toArray()[i].equals(item)) {
                queue.remove(queue.keySet().toArray()[i]);
                break;
            }
        }
    }
    
    public void DeleteItemForNumber(int number) {
        queue.remove(queue.keySet().toArray()[number]);
    }
    
    public Object getItemForKey(String Key) {
        return queue.get(Key);
    }
    
    public Object getItemForItem(Object item) {
        for (int i = 0; i < queue.values().toArray().length; i++) {
            if (queue.values().toArray()[i].equals(item)) {
                return item;
            }
        }
        return null;
    }
    
    public Object getItemForNumber(int number) {
        return queue.values().toArray()[number];
    }
    
    public String RandMd5Key() {
        byte[] bytesOfMessage = null, thedigest = null;
        String md5_key = "";
        
        while (true) {
            try  {
                int random_number = (int) (Math.random() * 1000000000); 
                bytesOfMessage = Integer.toString(random_number).getBytes("UTF-8");
                MessageDigest md = MessageDigest.getInstance("MD5");
                thedigest = md.digest(bytesOfMessage);
                for(int i = 0; i < thedigest.length;i++) {
                    md5_key = md5_key + (int)thedigest[i];
                }
            } catch (NoSuchAlgorithmException | UnsupportedEncodingException ex) {
                log_manager.CallError(ex.toString(), "QueueManager. Create md5_key", name_class);
            }
            boolean prov = true;
            for (String key : queue.keySet()) {
                if (key.equals(md5_key)) {
                    prov = false;
                    break;
                }
            }
            if (prov) 
                break;
        }
        return md5_key; 
    }
}
