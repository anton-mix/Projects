package TestSender;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import org.apache.commons.io.IOUtils;
import org.apache.http.HttpEntity;
import org.apache.http.HttpResponse;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpPost;
import org.apache.http.impl.client.HttpClients;

public class TestSender {

    static public String sendPost(String url, String[][] keys){
        HttpClient httpclient = HttpClients.createDefault();
        HttpResponse response;
        HttpEntity entity;
        InputStream instream = null;
        HttpPost httppost = new HttpPost(url);
        String result = "";
/*
        List<NameValuePair> params = new ArrayList<NameValuePair>(2);
        for (int i = 0; i < keys.length; i++){
            params.add(new BasicNameValuePair(keys[i][0], keys[i][1]));
        }
**/
        try {
            //httppost.setEntity(new UrlEncodedFormEntity(params, "UTF-8"));
            response = httpclient.execute(httppost);
            entity = response.getEntity();
            if (entity != null) {
                instream = entity.getContent();
                result = IOUtils.toString(instream, StandardCharsets.UTF_8);
            }
        } catch (Exception ex) {
            result = "{\"error\":{\"error_code\":99999,\"error_msg\":\"" + ex.toString() + "\",\"request_params\":\"\"}}";
        } finally {
            try {
                if (instream != null) {
                    instream.close();
                }
            } catch (IOException ex) {
                result = "{\"error\":{\"error_code\":99999,\"error_msg\":\"" + ex.toString() + "\",\"request_params\":\"\"}}";
            }
        }
        
        return result;
    }
}
