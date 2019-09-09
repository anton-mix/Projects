package services;


import ImageObjects.LoadImage;
import ImageObjects.LoadImageConverter;
import Managers.LogManager;
import Managers.MessageManager;
import MyEnum.ImageMicroserviseIdEnum;
import MyEnum.RestMailingEnum;
import MyEnum.RestShopEnum;
import MyEnum.RestUserEnum;
import MyInterface.RestInterface;
import RestObject.*;
import SQLClass.GlobalVariables;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.sun.jersey.multipart.FormDataParam;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.sql.SQLException;
import javax.ws.rs.Consumes;
import javax.ws.rs.PathParam;
import javax.ws.rs.GET;
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.Produces;
import javax.ws.rs.core.Context;
import javax.ws.rs.core.HttpHeaders;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;
import javax.ws.rs.core.UriInfo;


@Path("ManagerFiles")
public class ManagerFilesRest implements RestInterface {
    ResData answer;
    String name_class = "ManagerFilesRest";
    GlobalVariables global_variables;
    LogManager log_manager;
    MessageManager message_manager;
    Gson gson;
    
    public ManagerFilesRest() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerFilesRest");
        
        try {
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            builder.registerTypeAdapter(LoadImage.class, new LoadImageConverter());
            builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
            gson = builder.create();
            
            message_manager = new MessageManager(global_variables);
        }
        catch(Exception ex) {
            answer = log_manager.CallError(ex.toString(), "ManagerFilesRest. Init", name_class);
        } 
    }
    
    @GET
    @Path("/loadImage/{any: .*}")
    @Produces("image/jpg")
    public Response getImage(@Context UriInfo uriInfo) throws SQLException {
        File image;
        SQLResult RestRes;
        LoadImage BuffItem;
        ImageMicroserviseIdEnum imageEnum;
        int code_id = -1;
        try {
           BuffItem = gson.fromJson(uriInfo.getPath().split("/", 3)[2].replaceAll("/", "//"), LoadImage.class);
           imageEnum = ImageMicroserviseIdEnum.values()[BuffItem.getMicroservice_id()];
           switch (imageEnum) {
                case user_session_queue:
                    code_id = RestUserEnum.User_get_path_image.ordinal();
                    break;
                case shop_session_queue:
                    code_id = RestShopEnum.Product_get_path_image.ordinal();
                    break;
                case mailing_session_queue:
                    code_id = RestMailingEnum.Http_settings_get_path_image.ordinal();
                    break;
            }
           
            answer = CreateMessage(code_id, uriInfo.getPath().split("/", 3)[2].replaceAll("/", "//"), "getImage(UriInfo). Send message", ImageMicroserviseIdEnum.values()[BuffItem.getMicroservice_id()].name());
            if (!answer.getIs_success()) {
                throw new Exception("getImage(UriInfo). Send message");
            }
            RestRes = gson.fromJson(answer.getRes(), SQLResult.class);
            if (RestRes.getUser_id() == 1) {
                if(!"not Ok".equals(RestRes.getResult())) {
                    image = new File(RestRes.getResult());
                }
                else {
                    throw new Exception("Error");
                }
            } else {
                throw new Exception("Error");
            }
            return Response.ok(image, "image/jpg").header("Inline", "filename=\"" + image.getName() + "\"").build();
        } catch (Exception ex) {
            answer = log_manager.CallError(ex.toString(), "uploadImage(String). Exception", name_class);
            return Response.ok()
                .header("Access-Control-Allow-Origin", "*")
                .header("Access-Control-Allow-Methods", "GET, POST, DELETE, PUT")
                .entity(answer.getRes()).build();
        }
    }
    
    
    @GET
    @Path("/getCheckImage/{any: .*}")
    public Response getCheckImage(@Context UriInfo uriInfo) throws SQLException {
        LoadImage BuffItem;
        ImageMicroserviseIdEnum imageEnum;
        SQLResult res;
        int code_id = -1;
        try {
           BuffItem = gson.fromJson(uriInfo.getPath().split("/", 3)[2].replaceAll("/", "//"), LoadImage.class);
           imageEnum = ImageMicroserviseIdEnum.values()[BuffItem.getMicroservice_id()];
           switch (imageEnum) {
                case user_session_queue:
                    code_id = RestUserEnum.User_get_path_image.ordinal();
                    break;
                case shop_session_queue:
                    code_id = RestShopEnum.Product_get_path_image.ordinal();
                    break;
                case mailing_session_queue:
                    code_id = RestMailingEnum.Http_settings_get_path_image.ordinal();
                    break;
            }
           
            answer = CreateMessage(code_id, uriInfo.getPath().split("/", 3)[2].replaceAll("/", "//"), "getImage(UriInfo). Send message", ImageMicroserviseIdEnum.values()[BuffItem.getMicroservice_id()].name());
            if (!answer.getIs_success()) {
                throw new Exception("getImage(UriInfo). Send message");
            }
            res = gson.fromJson(answer.getRes(), SQLResult.class);
            if(res.getId() == 1) {
                answer.setRes(gson.toJson(new SQLResult("true"), SQLResult.class));
            }
            else if(res.getId() == 0) {
                answer.setRes(gson.toJson(new SQLResult("false"), SQLResult.class));
            }
        } catch (Exception ex) {
            answer = log_manager.CallError(ex.toString(), "uploadImage(String). Exception", name_class);
        }
        return Response.ok()
                .header("Access-Control-Allow-Origin", "*")
                .header("Access-Control-Allow-Methods", "GET, POST, DELETE, PUT")
                .entity(answer.getRes()).build();
    }
    

    @POST
    @Path("/uploadImage/{json}")
    @Consumes(MediaType.MULTIPART_FORM_DATA)
    public Response uploadImage(@Context HttpHeaders headers, @PathParam("json") String json, @FormDataParam("file") InputStream uploadedInputStream) {
        LoadImage BuffItem;
        SQLResult RestRes;
        ImageMicroserviseIdEnum imageEnum;
        int code_id = -1, size = 0;
        
        /*
        String res = "";
        for(int i = 0; i < headers.getRequestHeaders().size();i ++) {
            res = res + " name: " + headers.getRequestHeaders().keySet().toArray()[i] + " value: " + headers.getRequestHeaders().values().toArray()[i] + " <br>";
        }
        return res;
        */
        
        try {
            BuffItem = gson.fromJson(json, LoadImage.class);
            size = Integer.valueOf(headers.getHeaderString("Content-Length"));
            
            if (size > 5000000) {
                throw new Exception("Image to large");
            }
            
            imageEnum = ImageMicroserviseIdEnum.values()[BuffItem.getMicroservice_id()];
            switch (imageEnum) {
                case user_session_queue:
                    code_id = RestUserEnum.User_set_path_image.ordinal();
                    break;
                case shop_session_queue:
                    code_id = RestShopEnum.Product_set_path_image.ordinal();
                    break;
                case mailing_session_queue:
                    code_id = RestMailingEnum.Http_settings_set_path_image.ordinal();
                    break;
            }

            answer = CreateMessage(code_id, json, "uploadImage(String). Send message", ImageMicroserviseIdEnum.values()[BuffItem.getMicroservice_id()].name());
            if(!answer.getIs_success()) {
                throw new Exception("uploadImage(String). Send message");
            }
            
            RestRes = gson.fromJson(answer.getRes(), SQLResult.class);
            if (RestRes.getUser_id() != 1) {
                answer = new ResData(false, gson.toJson(new SQLResult("Error in load image"), SQLResult.class));
            } else {
                if(writeToFile(uploadedInputStream, RestRes.getResult(), Integer.valueOf(headers.getHeaderString("Content-Length")))) {
                    answer = new ResData(true, "{\"result\":\"Ok\"}");
                }
                else {
                    throw new Exception("Parse Error");
                }
            }
            
        } catch (Exception ex) {
            answer = log_manager.CallError(ex.toString(), "uploadImage(String). Exception", name_class);
        }
        
       return InitMessage(answer);
       //return res;
    }

    private Boolean writeToFile(InputStream uploadedInputStream, String uploadedFileLocation, int lenght) {
        byte[] image = new byte[lenght];
        int position = 0;
        int read = 0;
        byte[] bytes = new byte[1024];
        try {
            OutputStream out = new FileOutputStream(new File(uploadedFileLocation));
            while ((read = uploadedInputStream.read(bytes)) != -1) {
                System.arraycopy(bytes, 0, image, position, read);
                position = position + read;
            }
            read = position = 0;
            for (int i = 0; i < lenght; i++) {
                if (image[i] == '\n') {
                    read++;
                    if (read == 1) {
                        position = i;
                    }
                }
                if (read == 4) {
                    read = i;
                    break;
                }
            }
            out.write(image, read + 1, lenght - position - read - 6);
            out.flush();
            out.close();
        } catch (Exception e) {
            return false;
        }
        return true;
    }
    
    @Override
    public ResData CreateMessage(int code, String json, String ErrorMessage, String name_microservice) {
        try {
            JsonParser parser = new JsonParser();
            JsonObject json_obj = parser.parse(json).getAsJsonObject();
            String user_id = checkFildString(json_obj, "user_id");
            
            GsonBuilder builder = new GsonBuilder();
            builder.registerTypeAdapter(RestObject.class, new RestObjectConverter());
            return message_manager.call(builder.create().toJson(new RestObject(json_obj, user_id, code)), name_microservice);
        }
        catch (Exception ex) {
            return log_manager.CallError(ex.toString(), ErrorMessage, name_class);
        }
    }    
}
