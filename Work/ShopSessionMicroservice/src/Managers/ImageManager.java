package Managers;

import ImageObjects.LoadImage;
import ImageObjects.LoadImageConverter;
import MyEnum.ImageEnum;
import SQLClass.GlobalVariables;
import SQLClass.ResData;
import SQLClass.SQLResult;
import SQLClass.SQLResultConverter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import java.io.File;

public class ImageManager {
    LogManager log_manager;
    ProductsManager products_manager;
    GlobalVariables global_variables;
    Gson gson;
    
    String name_class = "LoadImageManager";
    
    public ImageManager(LogManager log_manager, GlobalVariables global_variables, ProductsManager products_manager) {
        this.log_manager = log_manager;
        this.global_variables = global_variables;
        this.products_manager = products_manager;
        
        GsonBuilder builder = new GsonBuilder();
        builder.registerTypeAdapter(LoadImage.class, new LoadImageConverter());
        builder.registerTypeAdapter(SQLResult.class, new SQLResultConverter());
        gson = builder.create();
    }
    
    public ResData setPathImage(LoadImage item) {
        ImageEnum imageEnum;
        ResData res;
        String BasePath = global_variables.getWindows_base_path_image();
        try {
            File file = new File(BasePath);

            if (!file.exists()) {
                BasePath = global_variables.getUnix_base_path_image();
            }
            
            imageEnum = ImageEnum.values()[item.getImage_id()];
            switch (imageEnum) {
                case logo:
                    if(item.getPart_path().size() != 1) {
                        throw new Exception("BadFormatException");
                    }
                    
                    BasePath = BasePath + item.getUser_id() + "//Products//" + item.getPart_path().get(0) + "//";
                    file = new File(BasePath);
                    if(!file.exists() && !file.mkdirs()) {
                        throw new Exception("FolderNotCreateException");
                    }
                    BasePath = BasePath + imageEnum.name() + ".jpg";
                    break;
            }
            
            res = new ResData(true, gson.toJson(new SQLResult(BasePath,1,1), SQLResult.class));
        }
        catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "getPathImage(LoadImage). Invalid Path", name_class);
        }
        return res;
    }
    
    public ResData getPathImage(LoadImage item) {
        ImageEnum imageEnum;
        ResData res;
        String BasePath = global_variables.getWindows_base_path_image(), CurPath = "";
        try {
            File file = new File(BasePath);

            if (!file.exists()) {
                BasePath = global_variables.getUnix_base_path_image();
            }
            
            imageEnum = ImageEnum.values()[item.getImage_id()];
            switch (imageEnum) {
                case logo:
                    if(item.getPart_path().size() != 1) {
                        throw new Exception("BadFormatException");
                    }
                    
                    CurPath = BasePath + item.getUser_id() + "//Products//" + item.getPart_path().get(0) + "//" + imageEnum.name() + ".jpg";
                    file = new File(CurPath);
                    if(!file.exists()) {
                       CurPath =  BasePath + "Default//ShopDefault//logo.jpg";
                    }
                    break;
            }
            
            res = new ResData(true, gson.toJson(new SQLResult(CurPath,1,1), SQLResult.class));
        }
        catch (Exception ex) {
            res = log_manager.CallError(ex.toString(), "getPathImage(LoadImage). Invalid Path", name_class);
        }
        return res;
    }
    
    public void setErrorHttp(String ErrorHttp) {
        this.log_manager.setErrorHttp(ErrorHttp);
    }
}
