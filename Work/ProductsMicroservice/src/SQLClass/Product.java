package SQLClass;

import java.sql.Timestamp;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

public class Product {
    transient int id = -1;
    transient String user_id = "";
    transient String product_type = "";
    transient String payment_type = "";
    transient String order_page = "";
    transient String product_name = "";
    transient String image = "";
    transient String product_description = "";
    transient String groups_after_prepayment = ""; 
    transient String message_topic_after_prepayment = "";
    transient String message_type_after_prepayment = ""; 
    transient String message_after_prepayment = "";
    transient String groups_after_payment = "";
    transient String thanks_page = "";
    transient String message_topic_after_payment = "";
    transient String message_type_after_payment = "";
    transient String message_after_payment = "";
    transient String getting_page = "";
    transient String pincodes = "";
    transient String vendor_code = "";
    transient String company = "";
    transient String country = "";
    transient String ordering = "";
    transient String bill_type = "";
    transient String employee_instruction = "";
    transient String product_type_for_add = "";
    transient String api_url_order = "";
    transient String api_url_prepayment = "";
    transient String api_url_payment = "";
    transient String api_url_cancel = "";
    transient String api_url_return = "";
    transient String api_url_alert_return = "";
    transient String personal_template = "";
    transient String some_payments = "";
    transient String show_in_catalogue = "";
    transient String tax = "";
    transient String show_for_partners = "";
    transient double value = -1.0;
    transient double prepayment_min = -1.0;
    transient double amount_of_expenses = -1.0;
    transient double percents_of_expenses = -1.0; 
    transient double comission = -1.0;
    transient double comission_percents = -1.0;
    transient double employee_reward = -1.0;
    transient double employee_reward_percents = -1.0;
    transient Timestamp bill_date;
    
    DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    Date date = dateFormat.parse("1970-01-02 00:00:00");
    long time = date.getTime();
    
    enum Product_enum {id, user_id, product_type,
      payment_type, order_page, product_name, image, product_description, groups_after_prepayment, message_topic_after_prepayment, 
      message_type_after_prepayment, message_after_prepayment, groups_after_payment, thanks_page, message_topic_after_payment,
      message_type_after_payment ,message_after_payment, getting_page, pincodes, vendor_code, company, country, ordering, bill_type,
      employee_instruction, product_type_for_add, api_url_order, api_url_prepayment, api_url_payment, api_url_cancel, api_url_return,
      api_url_alert_return, personal_template, some_payments, show_in_catalogue, tax, show_for_partners, value, prepayment_min,
     amount_of_expenses, percents_of_expenses, comission, comission_percents, employee_reward, employee_reward_percents, bill_date }
    
    public Product (int id, String user_id, String product_type, String payment_type, String order_page, String product_name,
            double value, String some_payments, double prepayment_min, String image, String show_in_catalogue, 
            String product_description, double amount_of_expenses, double percents_of_expenses, 
            String groups_after_prepayment, String message_topic_after_prepayment, 
            String message_type_after_prepayment, String message_after_prepayment,
            String groups_after_payment, String thanks_page, String message_topic_after_payment, 
            String message_type_after_payment, String message_after_payment, String getting_page, 
            String pincodes, String tax, String vendor_code, String company, String country,
            String ordering, String bill_type, Timestamp bill_date, double comission, 
            double comission_percents, String show_for_partners, double employee_reward, 
            double employee_reward_percents, String employee_instruction, String product_type_for_add, 
            String personal_template, String api_url_order, String api_url_prepayment, 
            String api_url_payment, String api_url_cancel, String api_url_return, String api_url_alert_return) throws ParseException  {
        this.id = id;
        this.user_id = user_id;
        this.product_type = product_type;
        this.payment_type = payment_type;
        this.order_page = order_page;
        this.product_name = product_name;
        this.value = value;
        this.some_payments = some_payments;
        this.prepayment_min = prepayment_min;
        this.image = image;
        this.show_in_catalogue = show_in_catalogue;
        this.product_description = product_description;
        this.amount_of_expenses = amount_of_expenses;
        this.percents_of_expenses = percents_of_expenses;
        this.groups_after_prepayment = groups_after_prepayment;
        this.message_topic_after_prepayment = message_topic_after_prepayment;
        this.message_type_after_prepayment = message_type_after_prepayment;
        this.message_after_prepayment = message_after_prepayment;
        this.groups_after_payment = groups_after_payment;
        this.thanks_page = thanks_page;
        this.message_topic_after_payment = message_topic_after_payment;
        this.message_type_after_payment = message_type_after_payment;
        this.message_after_payment = message_after_payment;
        this.getting_page = getting_page;
        this.pincodes = pincodes;
        this.tax = tax;
        this.vendor_code = vendor_code;
        this.company = company;
        this.country = country;
        this.ordering = ordering;
        this.bill_type = bill_type;
        if(bill_date == null) 
            this.bill_date = new Timestamp(time);
        else
            this.bill_date = bill_date;
        this.comission = comission;
        this.comission_percents = comission_percents;
        this.show_for_partners = show_for_partners;
        this.employee_reward = employee_reward;
        this.employee_reward_percents = employee_reward_percents;
        this.employee_instruction = employee_instruction;
        this.product_type_for_add = product_type_for_add;
        this.personal_template = personal_template;
        this.api_url_order = api_url_order;
        this.api_url_prepayment = api_url_prepayment;
        this.api_url_payment = api_url_payment;
        this.api_url_cancel = api_url_cancel;
        this.api_url_return = api_url_return;
        this.api_url_alert_return = api_url_alert_return;
    }
    
    public Product (String user_id, String product_type, String payment_type, String order_page, String product_name,
            double value, String some_payments, double prepayment_min, String image, String show_in_catalogue, 
            String product_description, double amount_of_expenses, double percents_of_expenses, 
            String groups_after_prepayment, String message_topic_after_prepayment, 
            String message_type_after_prepayment, String message_after_prepayment,
            String groups_after_payment, String thanks_page, String message_topic_after_payment, 
            String message_type_after_payment, String message_after_payment, String getting_page, 
            String pincodes, String tax, String vendor_code, String company, String country,
            String ordering, String bill_type, Timestamp bill_date, double comission, 
            double comission_percents, String show_for_partners, double employee_reward, 
            double employee_reward_percents, String employee_instruction, String product_type_for_add, 
            String personal_template, String api_url_order, String api_url_prepayment, 
            String api_url_payment, String api_url_cancel, String api_url_return, String api_url_alert_return) throws ParseException  {
        this.user_id = user_id;
        this.product_type = product_type;
        this.payment_type = payment_type;
        this.order_page = order_page;
        this.product_name = product_name;
        this.value = value;
        this.some_payments = some_payments;
        this.prepayment_min = prepayment_min;
        this.image = image;
        this.show_in_catalogue = show_in_catalogue;
        this.product_description = product_description;
        this.amount_of_expenses = amount_of_expenses;
        this.percents_of_expenses = percents_of_expenses;
        this.groups_after_prepayment = groups_after_prepayment;
        this.message_topic_after_prepayment = message_topic_after_prepayment;
        this.message_type_after_prepayment = message_type_after_prepayment;
        this.message_after_prepayment = message_after_prepayment;
        this.groups_after_payment = groups_after_payment;
        this.thanks_page = thanks_page;
        this.message_topic_after_payment = message_topic_after_payment;
        this.message_type_after_payment = message_type_after_payment;
        this.message_after_payment = message_after_payment;
        this.getting_page = getting_page;
        this.pincodes = pincodes;
        this.tax = tax;
        this.vendor_code = vendor_code;
        this.company = company;
        this.country = country;
        this.ordering = ordering;
        this.bill_type = bill_type;
        if(bill_date == null) 
            this.bill_date = new Timestamp(time);
        else
            this.bill_date = bill_date;
        this.comission = comission;
        this.comission_percents = comission_percents;
        this.show_for_partners = show_for_partners;
        this.employee_reward = employee_reward;
        this.employee_reward_percents = employee_reward_percents;
        this.employee_instruction = employee_instruction;
        this.product_type_for_add = product_type_for_add;
        this.personal_template = personal_template;
        this.api_url_order = api_url_order;
        this.api_url_prepayment = api_url_prepayment;
        this.api_url_payment = api_url_payment;
        this.api_url_cancel = api_url_cancel;
        this.api_url_return = api_url_return;
        this.api_url_alert_return = api_url_alert_return;
    }
    
    
    public Product(String json) throws ParseException  {
        DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        this.bill_date = new Timestamp(time);
        List<String[]> jsonString = GlobalVariables.Parse_json(json);

        for(int i = 0; i < jsonString.size(); i ++) {

            Product_enum product_enum = Product_enum.valueOf(jsonString.get(i)[0]);
            
            switch(product_enum) {
                case id:
                    this.id = Integer.valueOf(jsonString.get(i)[1]);
                    break;
                case user_id:
                    this.user_id = jsonString.get(i)[1];
                    break;
                case product_type:
                    this.product_type = jsonString.get(i)[1];
                    break;
                case payment_type:
                    this.payment_type = jsonString.get(i)[1];
                    break;
                case order_page:
                    this.order_page = jsonString.get(i)[1];
                    break;
                case product_name:
                    this.product_name = jsonString.get(i)[1];
                    break;
                case value:
                    this.value = Double.valueOf(jsonString.get(i)[1]);
                    break;
                case some_payments:
                    this.some_payments = jsonString.get(i)[1];
                    break;
                case prepayment_min:
                    this.prepayment_min = Double.valueOf(jsonString.get(i)[1]);
                    break;
                case image:
                    this.image = jsonString.get(i)[1];
                    break;
                case show_in_catalogue:
                    this.show_in_catalogue = jsonString.get(i)[1];
                    break;
                case product_description:
                    this.product_description = jsonString.get(i)[1];
                    break;
                case amount_of_expenses:
                    this.amount_of_expenses = Double.valueOf(jsonString.get(i)[1]);
                    break;
                case percents_of_expenses:
                    this.percents_of_expenses = Double.valueOf(jsonString.get(i)[1]);
                    break;
                case groups_after_prepayment:
                    this.groups_after_prepayment = jsonString.get(i)[1];
                    break;
                case message_topic_after_prepayment:
                    this.message_topic_after_prepayment = jsonString.get(i)[1];
                    break;
                case message_type_after_prepayment:
                    this.message_type_after_prepayment = jsonString.get(i)[1];
                    break;
                case message_after_prepayment:
                    this.message_after_prepayment = jsonString.get(i)[1];
                    break;
                case groups_after_payment:
                    this.groups_after_payment = jsonString.get(i)[1];
                    break;
                case thanks_page:
                    this.thanks_page = jsonString.get(i)[1];
                    break;
                case message_topic_after_payment:
                    this.message_topic_after_payment = jsonString.get(i)[1];
                    break;
                case message_type_after_payment:
                    this.message_type_after_payment = jsonString.get(i)[1];
                    break;
                case message_after_payment:
                    this.message_after_payment = jsonString.get(i)[1];
                    break;
                case getting_page:
                    this.getting_page = jsonString.get(i)[1];
                    break;
                case pincodes:
                    this.pincodes = jsonString.get(i)[1];
                    break;
                case tax:
                    this.tax = jsonString.get(i)[1];
                    break;
                case vendor_code:
                    this.vendor_code = jsonString.get(i)[1];
                    break;
                case company:
                    this.company = jsonString.get(i)[1];
                    break;
                case country:
                    this.country = jsonString.get(i)[1];
                    break;
                case ordering:
                    this.ordering = jsonString.get(i)[1];
                    break;
                case bill_type:
                    this.bill_type = jsonString.get(i)[1];
                    break;
                case bill_date:
                    try {
                    Date date = dateFormat.parse(jsonString.get(i)[1]);
                    long time = date.getTime();
                    this.bill_date = new Timestamp(time);
                    }
                    catch(ParseException ex) {
                        System.out.println("Time error" + ex.toString());
                        throw ex;
                    }
                    break;
                case comission:
                    this.comission = Double.valueOf(jsonString.get(i)[1]);
                    break;
                case comission_percents:
                    this.comission_percents = Double.valueOf(jsonString.get(i)[1]);
                    break;
                case show_for_partners:
                    this.show_for_partners = jsonString.get(i)[1];
                    break;
                case employee_reward:
                    this.employee_reward = Double.valueOf(jsonString.get(i)[1]);
                    break;
                case employee_reward_percents:
                    this.employee_reward_percents = Double.valueOf(jsonString.get(i)[1]);
                    break;
                case employee_instruction:
                    this.employee_instruction = jsonString.get(i)[1];
                    break;
                case product_type_for_add:
                    this.product_type_for_add = jsonString.get(i)[1];
                    break;
                case personal_template:
                    this.personal_template = jsonString.get(i)[1];
                    break;
                case api_url_order:
                    this.api_url_order = jsonString.get(i)[1];
                    break;
                case api_url_prepayment:
                    this.api_url_prepayment = jsonString.get(i)[1];
                    break;
                case api_url_payment:
                    this.api_url_payment = jsonString.get(i)[1];
                    break;
                case api_url_cancel:
                    this.api_url_cancel = jsonString.get(i)[1];
                    break;
                case api_url_return:
                    this.api_url_return = jsonString.get(i)[1];
                    break;
                case api_url_alert_return:
                    this.api_url_alert_return = jsonString.get(i)[1];
                    break;
            }
        }
    }
    
    public String[] getArray() {
        return new String[] {String.valueOf(id), user_id, product_type, payment_type, order_page, product_name, String.valueOf(value), 
            some_payments, String.valueOf(prepayment_min), image, show_in_catalogue, product_description, 
            String.valueOf(amount_of_expenses), String.valueOf(percents_of_expenses), groups_after_prepayment, message_topic_after_prepayment, 
            message_type_after_prepayment, message_after_prepayment, groups_after_payment, thanks_page, message_topic_after_payment, 
            message_type_after_payment, message_after_payment, getting_page, pincodes, tax, vendor_code, company, country,
            ordering, bill_type, bill_date.toString(), String.valueOf(comission), String.valueOf(comission_percents), 
            show_for_partners, String.valueOf(employee_reward), String.valueOf(employee_reward_percents), 
            employee_instruction, product_type_for_add, personal_template, api_url_order, api_url_prepayment, api_url_payment, 
            api_url_cancel, api_url_return, api_url_alert_return};
    }

    public String[] getSqlArray() {
        return new String[] {String.valueOf(id),
            "'" + user_id.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + product_type.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + payment_type.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + order_page.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + product_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            String.valueOf(value), 
            "'" + some_payments.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"")  + "'", 
            String.valueOf(prepayment_min),
            "'" + image.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + show_in_catalogue.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",        
            "'" + product_description.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            String.valueOf(amount_of_expenses), 
            String.valueOf(percents_of_expenses),
            "'" + groups_after_prepayment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + message_topic_after_prepayment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + message_type_after_prepayment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + message_after_prepayment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + groups_after_payment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + thanks_page.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + message_topic_after_payment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + message_type_after_payment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + message_after_payment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + getting_page.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + pincodes.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + tax.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + vendor_code.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + company.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + country.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + ordering.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + bill_type.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + bill_date.toString() + "'",
            String.valueOf(comission), 
            String.valueOf(comission_percents), 
            "'" + show_for_partners.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            String.valueOf(employee_reward), 
            String.valueOf(employee_reward_percents),
            "'" + employee_instruction.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + product_type_for_add.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + personal_template.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + api_url_order.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + api_url_prepayment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + api_url_payment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + api_url_cancel.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + api_url_return.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + api_url_alert_return.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"};
    }
    
    public String getSqlStringArray() {
        return "null ," +
            "'" + user_id.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + product_type.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + payment_type.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + order_page.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + product_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            String.valueOf(value) + "," +
            "'" + some_payments.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            String.valueOf(prepayment_min) + "," +
            "'" + image.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + show_in_catalogue.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +       
            "'" + product_description.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            String.valueOf(amount_of_expenses) + "," +
            String.valueOf(percents_of_expenses) + "," +
            "'" + groups_after_prepayment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + message_topic_after_prepayment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + message_type_after_prepayment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," + 
            "'" + message_after_prepayment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + groups_after_payment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + thanks_page.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + message_topic_after_payment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + message_type_after_payment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," + 
            "'" + message_after_payment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," + 
            "'" + getting_page.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + pincodes.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + tax.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + vendor_code.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + company.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + country.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + ordering.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," + 
            "'" + bill_type.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + bill_date.toString() + "'," +
            String.valueOf(comission) + "," +
            String.valueOf(comission_percents) + "," +
            "'" + show_for_partners.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            String.valueOf(employee_reward) + "," +
            String.valueOf(employee_reward_percents) + "," +
            "'" + employee_instruction.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + product_type_for_add.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," + 
            "'" + personal_template.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + api_url_order.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + api_url_prepayment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + api_url_payment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + api_url_cancel.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + api_url_return.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + api_url_alert_return.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'";
    }
        
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "product_type", "payment_type", "order_page", "product_name", "value", "some_payments", "prepayment_min", 
            "image", "show_in_catalogue", "product_description", "amount_of_expenses", "percents_of_expenses", "groups_after_prepayment",
            "message_topic_after_prepayment", "message_type_after_prepayment", "message_after_prepayment", "groups_after_payment", "thanks_page",
            "message_topic_after_payment", "message_type_after_payment", "message_after_payment", "getting_page", "pincodes", "tax", "vendor_code", "company", "country",
            "ordering", "bill_type", "bill_date", "comission", "comission_percents", "show_for_partners",  "employee_reward", "employee_reward_percents",
            "employee_instruction", "product_type_for_add", "personal_template", "api_url_order", "api_url_prepayment", "api_url_payment", "api_url_cancel",
            "api_url_return", "api_url_alert_return"};
    }
    
    public String getKeyString() {
        return "id, user_id, product_type, payment_type, order_page, product_name, value, some_payments, prepayment_min, " + 
            "image, show_in_catalogue, product_description, amount_of_expenses, percents_of_expenses, groups_after_prepayment, " +
            "message_topic_after_prepayment, message_type_after_prepayment, message_after_prepayment, groups_after_payment, thanks_page, " +
            "message_topic_after_payment, message_type_after_payment, message_after_payment, getting_page, pincodes, tax, vendor_code, company, country, " +
            "ordering, bill_type, bill_date, comission, comission_percents, show_for_partners,  employee_reward, employee_reward_percents, " +
            "employee_instruction, product_type_for_add, personal_template, api_url_order, api_url_prepayment, api_url_payment, api_url_cancel, "+
            "api_url_return, api_url_alert_return"; 
    }
        
    public String getJson() {
        String[] buffItemString = getArray();
        String[] buffItemKeyString = getKeyArray();
        String Result = "{";
        
        for(int i = 0;i < buffItemString.length;i++) {
            Result = Result +  "\"" + buffItemKeyString[i] + "\":\"" + buffItemString[i] + "\"";
            if(i != buffItemString.length-1)
                    Result = Result + ",";
            }
        return Result +  "}";
    }
    
    public int getId() {
        return id;
    }

    public String getProduct_type() {
        return product_type;
    }

    public String getPayment_type() {
        return payment_type;
    }

    public String getOrder_page() {
        return order_page;
    }

    public String getProduct_name() {
        return product_name;
    }

    public String getImage() {
        return image;
    }

    public String getProduct_description() {
        return product_description;
    }

    public String getGroups_after_prepayment() {
        return groups_after_prepayment;
    }

    public String getMessage_topic_after_prepayment() {
        return message_topic_after_prepayment;
    }

    public String getMessage_type_after_prepayment() {
        return message_type_after_prepayment;
    }

    public String getMessage_after_prepayment() {
        return message_after_prepayment;
    }

    public String getGroups_after_payment() {
        return groups_after_payment;
    }

    public String getThanks_page() {
        return thanks_page;
    }

    public String getMessage_topic_after_payment() {
        return message_topic_after_payment;
    }

    public String getMessage_type_after_payment() {
        return message_type_after_payment;
    }

    public String getMessage_after_payment() {
        return message_after_payment;
    }

    public String getGetting_page() {
        return getting_page;
    }

    public String getPincodes() {
        return pincodes;
    }

    public String getVendor_code() {
        return vendor_code;
    }

    public String getCompany() {
        return company;
    }

    public String getCountry() {
        return country;
    }

    public String getOrdering() {
        return ordering;
    }

    public String getBill_type() {
        return bill_type;
    }

    public String getEmployee_instruction() {
        return employee_instruction;
    }

    public String getProduct_type_for_add() {
        return product_type_for_add;
    }

    public String getApi_url_order() {
        return api_url_order;
    }

    public String getApi_url_prepayment() {
        return api_url_prepayment;
    }

    public String getApi_url_payment() {
        return api_url_payment;
    }

    public String getApi_url_cancel() {
        return api_url_cancel;
    }

    public String getApi_url_return() {
        return api_url_return;
    }

    public String getApi_url_alert_return() {
        return api_url_alert_return;
    }

    public double getValue() {
        return value;
    }

    public double getPrepayment_min() {
        return prepayment_min;
    }

    public double getAmount_of_expenses() {
        return amount_of_expenses;
    }

    public double getPercents_of_expenses() {
        return percents_of_expenses;
    }

    public double getComission() {
        return comission;
    }

    public double getComission_percents() {
        return comission_percents;
    }

    public double getEmployee_reward() {
        return employee_reward;
    }

    public double getEmployee_reward_percents() {
        return employee_reward_percents;
    }

    public Timestamp getBill_date() {
        return bill_date;
    }

    public String isPersonal_template() {
        return personal_template;
    }

    public String isSome_payments() {
        return some_payments;
    }

    public String isShow_in_catalogue() {
        return show_in_catalogue;
    }

    public String isTax() {
        return tax;
    }

    public String isShow_for_partners() {
        return show_for_partners;
    }

    public void setId(int id) {
        this.id = id;
    }

    public void setProduct_type(String product_type) {
        this.product_type = product_type;
    }

    public void setPayment_type(String payment_type) {
        this.payment_type = payment_type;
    }

    public void setOrder_page(String order_page) {
        this.order_page = order_page;
    }

    public void setProduct_name(String product_name) {
        this.product_name = product_name;
    }

    public void setImage(String image) {
        this.image = image;
    }

    public void setProduct_description(String product_description) {
        this.product_description = product_description;
    }

    public void setGroups_after_prepayment(String groups_after_prepayment) {
        this.groups_after_prepayment = groups_after_prepayment;
    }

    public void setMessage_topic_after_prepayment(String message_topic_after_prepayment) {
        this.message_topic_after_prepayment = message_topic_after_prepayment;
    }

    public void setMessage_type_after_prepayment(String message_type_after_prepayment) {
        this.message_type_after_prepayment = message_type_after_prepayment;
    }

    public void setMessage_after_prepayment(String message_after_prepayment) {
        this.message_after_prepayment = message_after_prepayment;
    }

    public void setGroups_after_payment(String groups_after_payment) {
        this.groups_after_payment = groups_after_payment;
    }

    public void setThanks_page(String thanks_page) {
        this.thanks_page = thanks_page;
    }

    public void setMessage_topic_after_payment(String message_topic_after_payment) {
        this.message_topic_after_payment = message_topic_after_payment;
    }

    public void setMessage_type_after_payment(String message_type_after_payment) {
        this.message_type_after_payment = message_type_after_payment;
    }

    public void setMessage_after_payment(String message_after_payment) {
        this.message_after_payment = message_after_payment;
    }

    public void setGetting_page(String getting_page) {
        this.getting_page = getting_page;
    }

    public void setPincodes(String pincodes) {
        this.pincodes = pincodes;
    }

    public void setVendor_code(String vendor_code) {
        this.vendor_code = vendor_code;
    }

    public void setCompany(String company) {
        this.company = company;
    }

    public void setCountry(String country) {
        this.country = country;
    }

    public void setOrdering(String ordering) {
        this.ordering = ordering;
    }

    public void setBill_type(String bill_type) {
        this.bill_type = bill_type;
    }

    public void setEmployee_instruction(String employee_instruction) {
        this.employee_instruction = employee_instruction;
    }

    public void setProduct_type_for_add(String product_type_for_add) {
        this.product_type_for_add = product_type_for_add;
    }

    public void setApi_url_order(String api_url_order) {
        this.api_url_order = api_url_order;
    }

    public void setApi_url_prepayment(String api_url_prepayment) {
        this.api_url_prepayment = api_url_prepayment;
    }

    public void setApi_url_payment(String api_url_payment) {
        this.api_url_payment = api_url_payment;
    }

    public void setApi_url_cancel(String api_url_cancel) {
        this.api_url_cancel = api_url_cancel;
    }

    public void setApi_url_return(String api_url_return) {
        this.api_url_return = api_url_return;
    }

    public void setApi_url_alert_return(String api_url_alert_return) {
        this.api_url_alert_return = api_url_alert_return;
    }

    public void setValue(double value) {
        this.value = value;
    }

    public void setPrepayment_min(double prepayment_min) {
        this.prepayment_min = prepayment_min;
    }

    public void setAmount_of_expenses(double amount_of_expenses) {
        this.amount_of_expenses = amount_of_expenses;
    }

    public void setPercents_of_expenses(double percents_of_expenses) {
        this.percents_of_expenses = percents_of_expenses;
    }

    public void setComission(double comission) {
        this.comission = comission;
    }

    public void setComission_percents(double comission_percents) {
        this.comission_percents = comission_percents;
    }

    public void setEmployee_reward(double employee_reward) {
        this.employee_reward = employee_reward;
    }

    public void setEmployee_reward_percents(double employee_reward_percents) {
        this.employee_reward_percents = employee_reward_percents;
    }

    public void setBill_date(Timestamp bill_date) {
        this.bill_date = bill_date;
    }

    public void setPersonal_template(String personal_template) {
        this.personal_template = personal_template;
    }

    public void setSome_payments(String some_payments) {
        this.some_payments = some_payments;
    }

    public void setShow_in_catalogue(String show_in_catalogue) {
        this.show_in_catalogue = show_in_catalogue;
    }

    public void setTax(String tax) {
        this.tax = tax;
    }

    public void setShow_for_partners(String show_for_partners) {
        this.show_for_partners = show_for_partners;
    }
    
    public String getUser_id() {
        return user_id;
    }

    public void setUser_id(String user_id) {
        this.user_id = user_id;
    }
}