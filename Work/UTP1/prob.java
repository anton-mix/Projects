/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package services;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.sql.SQLException;
import java.util.Base64;
import javax.ws.rs.Consumes;
import javax.ws.rs.FormParam;
import javax.ws.rs.PathParam;
import javax.ws.rs.Produces;
import javax.ws.rs.GET;
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;


@Path("generic")
public class prob {
/*
    @GET
    @Path("/img3/{id}")
    @Produces("image/jpg")
    public Response getFile(@PathParam("id") String id) throws SQLException {
        File file = new File("C://image//" + id + ".jpg");
        return Response.ok(file, "image/jpg").header("Inline", "filename=\"" + file.getName() + "\"")
                .build();
    }
    @POST
    @Path("/upload/{primaryKey}")
    @Consumes(MediaType.APPLICATION_FORM_URLENCODED)
    @Produces("image/jpg")
    public String uploadImage(@FormParam("image") String image, @PathParam("primaryKey") String primaryKey) throws SQLException, FileNotFoundException {
        String result = "false";
        FileOutputStream fos;

        fos = new FileOutputStream("C://" + primaryKey + ".jpg");

        // decode Base64 String to image
        try {

            byte byteArray[] = Base64.getMimeDecoder().decode(image);
            fos.write(byteArray);

            result = "true";
            fos.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
*/
}