/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package services;

import MyEnum.PaymentsOperations;
import Managers.LogManager;
import Managers.MessageManager;
import MyInterface.RestInterface;
import RestObject.RestObject;
import RestObject.RestObjectConverter;
import SQLClass.GlobalVariables;
import SQLClass.ResData;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import java.io.IOException;
import java.util.concurrent.TimeoutException;
import javax.ws.rs.Consumes;
import javax.ws.rs.FormParam;
import javax.ws.rs.GET;
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.core.Response;
import javax.servlet.http.HttpServletRequest;
import javax.ws.rs.core.Context;
import javax.ws.rs.core.MultivaluedMap;
import java.util.List;
import java.util.HashMap;
import java.util.Map;

/**
 * This class represents REST communication with PaymentsMicroservice
 * @author Orachigami
 */
@Path("ManagerPayments")
public class ManagerPayments implements RestInterface {
    ResData answer;
    String name_class = "ManagerPayments";
    GlobalVariables global_variables;
    MessageManager message_manager;
    LogManager log_manager;

    public ManagerPayments() {
        global_variables = new GlobalVariables();
        log_manager = new LogManager(global_variables, "ManagerPayments");
        try {
            message_manager = new MessageManager(global_variables);
        } catch (IOException | TimeoutException e) {
            log_manager.CallError(e.toString(), "ManagerPayments(). Init", name_class);
        }
    }

    private static final String[] IP_HEADER_CANDIDATES = {
        "X-Forwarded-For",
        "Proxy-Client-IP",
        "WL-Proxy-Client-IP",
        "HTTP_X_FORWARDED_FOR",
        "HTTP_X_FORWARDED",
        "HTTP_X_CLUSTER_CLIENT_IP",
        "HTTP_CLIENT_IP",
        "HTTP_FORWARDED_FOR",
        "HTTP_FORWARDED",
        "HTTP_VIA",
        "REMOTE_ADDR"
    };

    @GET
    @Path("getIp")
    public Response getIp(@Context HttpServletRequest request) {
        for (String header : IP_HEADER_CANDIDATES) {
            String ip = request.getHeader(header);
            if (ip != null && ip.length() != 0 && !"unknown".equalsIgnoreCase(ip)) {
                return InitMessage(new ResData(true, "Your IP: " + ip));
            }
        }
        return InitMessage(new ResData(true, "Your IP: " + request.getRemoteAddr()));
    }

    @POST
    @Path("notifyYandex")
    public Response notifyYandex(MultivaluedMap<String, String> formParams) {
        Map<String, String> data = new HashMap(22);
        for (Map.Entry<String, List<String>> param : formParams.entrySet()) {
            data.put(param.getKey(), param.getValue().get(0));
        }
        answer = CreateMessage(PaymentsOperations.NOTIFY_YANDEX, new Gson().toJson(data), "notifyYandex(String). Send message", "payments_queue");
        return InitMessage(answer);
    }

    @GET
    @Path("notifyPayu")
    public Response notifyPayu() {
        return Response.ok().build();
    }

    @POST
    @Path("notifyPayu")
    public Response notifyPayu(MultivaluedMap<String, String> formParams) {
        Gson gson = new Gson();
        Map<String, String> params = new HashMap(32);
        Map<String, List<String>> array_params = new HashMap(22);
        for (Map.Entry<String, List<String>> param : formParams.entrySet()) {
            if (param.getKey().endsWith("[]")) {
                params.put(param.getKey(), param.getValue().get(0));
            } else {
                array_params.put(param.getKey(), param.getValue());
            }
        }
        String part1 = gson.toJson(params);
        String part2 = gson.toJson(array_params);
        String json = part1.substring(0, part1.length() - 1) + "," + part2.substring(1);

        answer = CreateMessage(PaymentsOperations.NOTIFY_PAYU, json, "notifyPayu(String). Send message", "payments_queue");
        return InitMessage(answer);
    }

    @POST
    @Consumes("application/x-www-form-urlencoded")
    @Path("order")
    public Response order(@FormParam("json") String json) {
        try {
            answer = CreateMessage(PaymentsOperations.ORDER, json, "order (String). Send message", "payments_queue");
        } catch (Exception e) {
            return Response.ok(e.toString()).build();
        }
        return InitMessage(answer);
    }

    @POST
    @Path("addPayment")
    @Consumes("application/x-www-form-urlencoded")
    public Response addPayment(@FormParam("json") String json) {
        answer = CreateMessage(PaymentsOperations.ADD_PAYMENT, json, "addPayment(String). Send message", "payments_queue");
        return InitMessage(answer);
    }

    @POST
    @Path("getPayments")
    @Consumes("application/x-www-form-urlencoded")
    public Response getPayments(@FormParam("json") String json) {
        answer = CreateMessage(PaymentsOperations.GET_PAYMENTS, json, "getPayments(String). Send message", "payments_queue");
        return InitMessage(answer);
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
        } catch (Exception ex) {
            return log_manager.CallError(ex.toString(), ErrorMessage, name_class);
        }
    }
}
