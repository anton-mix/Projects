package Managers;

import SQLClass.GlobalVariables;
import java.io.IOException;
import net.spy.memcached.AddrUtil;
import net.spy.memcached.ConnectionFactoryBuilder;
import net.spy.memcached.MemcachedClient;
import net.spy.memcached.auth.AuthDescriptor;
import net.spy.memcached.auth.PlainCallbackHandler;

public class MemCachedManager {
    MemcachedClient mc;
    LogManager log_manager;
    GlobalVariables global_variables;
    String name_class = "MemCachedManager";
    
    public MemCachedManager(LogManager log_manager, GlobalVariables global_variables) {
        AuthDescriptor ad = new AuthDescriptor(new String[]{"PLAIN"}, new PlainCallbackHandler(global_variables.getDb_mem_user_name(), global_variables.getDb_mem_password()));
        try {
            mc = new MemcachedClient(
                    new ConnectionFactoryBuilder()
                            .setProtocol(ConnectionFactoryBuilder.Protocol.BINARY)
                            .setAuthDescriptor(ad).build(),
                    AddrUtil.getAddresses(global_variables.getDb_ip_mem() + ":" + global_variables.getDb_port_mem()));
        } catch (Exception ex) {
            log_manager.CallError(ex.toString(), "MemCachedManager. Error on init", name_class);
            System.out.println(ex.toString());
        }
    }
    
    public void SetItem(String key, String value) {
        mc.set(key, 0, value);
    }
    
    public String GetItem(String key) {
        return String.valueOf(mc.get(key));
    }
}
