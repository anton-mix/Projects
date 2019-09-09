package SQLClass.Product;

import MyInterface.ModelInterface;
import java.sql.Timestamp;
import java.text.ParseException;

public class ProductModel implements ModelInterface {
    int id = -1;
    int user_id = -1;
    int group_id = -1;
    String product_type = "";
    int payment_type = -1;
    String order_page = "";
    String product_name = "";
    String identifier = "";
    String image = "";
    String product_description = "";
    String groups_after_prepayment = ""; 
    String message_topic_after_prepayment = "";
    String message_type_after_prepayment = ""; 
    String message_after_prepayment = "";
    String groups_after_payment = "";
    String thanks_page = "";
    String message_topic_after_payment = "";
    String message_type_after_payment = "";
    String message_after_payment = "";
    String getting_page = "";
    String pincodes = "";
    String vendor_code = "";
    String company = "";
    String country = "";
    String ordering = "";
    String bill_type = "";
    String employee_instruction = "";
    String personal_template = "";
    String some_payments = "";
    String show_in_catalogue = "";
    String tax = "";
    String show_for_partners = "";
    double value = -1.0;
    double prepayment_min = -1.0;
    double amount_of_expenses = -1.0;
    double percents_of_expenses = -1.0; 
    double comission = -1.0;
    double comission_percents = -1.0;
    double employee_reward = -1.0;
    double employee_reward_percents = -1.0;
    String personal_topic = ""; 
    String personal_text = "";
    Timestamp bill_date;
    int reminder_payment_id = -1;
    int reminder_surcharge_id = -1;
    
    public ProductModel() {
        this.id = -1;
        this.user_id = -1;
        this.group_id = -1;
        this.product_type = "";
        this.payment_type = -1;
        this.order_page = "";
        this.product_name = "";
        this.identifier = "";
        this.image = "";
        this.product_description = "";
        this.groups_after_prepayment = "";
        this.message_topic_after_prepayment = "";
        this.message_type_after_prepayment = "";
        this.message_after_prepayment = "";
        this.groups_after_payment = "";
        this.thanks_page = "";
        this.message_topic_after_payment = "";
        this.message_type_after_payment = "";
        this.message_after_payment = "";
        this.getting_page = "";
        this.pincodes = "";
        this.vendor_code = "";
        this.company = "";
        this.country = "";
        this.ordering = "";
        this.bill_type = "";
        this.employee_instruction = "";
        this.personal_template = "";
        this.some_payments = "";
        this.show_in_catalogue = "";
        this.tax = "";
        this.show_for_partners = "";
        this.value = -1.0;
        this.prepayment_min = -1.0;
        this.amount_of_expenses = -1.0;
        this.percents_of_expenses = -1.0;
        this.comission = -1.0;
        this.comission_percents = -1.0;
        this.employee_reward = -1.0;
        this.employee_reward_percents = -1.0;
        this.personal_topic = "";
        this.personal_text = "";
        this.bill_date = Timestamp.valueOf("1970-01-02 0:0:0");
        this.reminder_payment_id = -1;
        this.reminder_surcharge_id = -1;
    }
    
    
    public ProductModel (int id, int user_id, int group_id, String product_type, int payment_type, String order_page, String product_name, 
            String identifier, double value, String some_payments, double prepayment_min, String image, String show_in_catalogue, 
            String product_description, double amount_of_expenses, double percents_of_expenses, String groups_after_prepayment, 
            String message_topic_after_prepayment, String message_type_after_prepayment, String message_after_prepayment,
            String groups_after_payment, String thanks_page, String message_topic_after_payment, String message_type_after_payment, 
            String message_after_payment, String getting_page, String pincodes, String tax, String vendor_code, String company, String country,
            String ordering, String bill_type, Timestamp bill_date, double comission, double comission_percents, String show_for_partners, 
            double employee_reward, double employee_reward_percents, String employee_instruction,String personal_template, 
            String personal_topic, String personal_text, int reminder_payment_id, int reminder_surcharge_id) {
        this.id = id;
        this.user_id = user_id;
        this.group_id = group_id;
        this.product_type = product_type;
        this.payment_type = payment_type;
        this.order_page = order_page;
        this.product_name = product_name;
        this.identifier = identifier;
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
        this.bill_date = bill_date;
        this.comission = comission;
        this.comission_percents = comission_percents;
        this.show_for_partners = show_for_partners;
        this.employee_reward = employee_reward;
        this.employee_reward_percents = employee_reward_percents;
        this.employee_instruction = employee_instruction;
        this.personal_template = personal_template;
        this.personal_topic = personal_topic;
        this.reminder_payment_id = reminder_payment_id;
        this.reminder_surcharge_id = reminder_surcharge_id;
    }
    
    public ProductModel (int user_id, String product_type, int payment_type, String order_page, String product_name, String identifier,
            double value, String some_payments, double prepayment_min, String image, String show_in_catalogue, String product_description, 
            double amount_of_expenses, double percents_of_expenses, String groups_after_prepayment, String message_topic_after_prepayment, 
            String message_type_after_prepayment, String message_after_prepayment, String groups_after_payment, String thanks_page, 
            String message_topic_after_payment, String message_type_after_payment, String message_after_payment, String getting_page, 
            String pincodes, String tax, String vendor_code, String company, String country, String ordering, String bill_type, 
            Timestamp bill_date, double comission, double comission_percents, String show_for_partners, double employee_reward, 
            double employee_reward_percents, String employee_instruction, String personal_template, String personal_topic, 
            String personal_text, int reminder_payment_id, int reminder_surcharge_id) throws ParseException  {
        this.user_id = user_id;
        this.product_type = product_type;
        this.payment_type = payment_type;
        this.order_page = order_page;
        this.product_name = product_name;
        this.identifier = identifier;
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
        this.bill_date = bill_date;
        this.comission = comission;
        this.comission_percents = comission_percents;
        this.show_for_partners = show_for_partners;
        this.employee_reward = employee_reward;
        this.employee_reward_percents = employee_reward_percents;
        this.employee_instruction = employee_instruction;
        this.personal_template = personal_template;
        this.personal_text = personal_text;
        this.personal_topic = personal_topic;
        this.reminder_payment_id = reminder_payment_id;
        this.reminder_surcharge_id = reminder_surcharge_id;
    }
    
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), String.valueOf(group_id),product_type, String.valueOf(payment_type), order_page, product_name, identifier, String.valueOf(value), 
            some_payments, String.valueOf(prepayment_min), image, show_in_catalogue, product_description, 
            String.valueOf(amount_of_expenses), String.valueOf(percents_of_expenses), groups_after_prepayment, message_topic_after_prepayment, 
            message_type_after_prepayment, message_after_prepayment, groups_after_payment, thanks_page, message_topic_after_payment, 
            message_type_after_payment, message_after_payment, getting_page, pincodes, tax, vendor_code, company, country, ordering, bill_type, 
            bill_date.toString(), String.valueOf(comission), String.valueOf(comission_percents), show_for_partners, 
            String.valueOf(employee_reward), String.valueOf(employee_reward_percents), employee_instruction, personal_template, personal_text, 
            personal_topic, String.valueOf(reminder_payment_id), String.valueOf(reminder_surcharge_id)};
    }

    @Override
    public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id),
            String.valueOf(user_id),
            String.valueOf(group_id),
            "'" + product_type.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            String.valueOf(payment_type),
            "'" + order_page.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + product_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + identifier.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
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
            "'" + personal_template.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + personal_topic.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + personal_text.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            String.valueOf(reminder_payment_id), 
            String.valueOf(reminder_surcharge_id)
        };
    }
        
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "group_id", "product_type", "payment_type", "order_page", "product_name", "identifier", "value", "some_payments", "prepayment_min", 
            "image", "show_in_catalogue", "product_description", "amount_of_expenses", "percents_of_expenses", "groups_after_prepayment",
            "message_topic_after_prepayment", "message_type_after_prepayment", "message_after_prepayment", "groups_after_payment", "thanks_page",
            "message_topic_after_payment", "message_type_after_payment", "message_after_payment", "getting_page", "pincodes", "tax", "vendor_code", "company", "country",
            "ordering", "bill_type", "bill_date", "comission", "comission_percents", "show_for_partners",  "employee_reward", "employee_reward_percents",
            "employee_instruction", "personal_template", "personal_topic", "personal_text", "reminder_payment_id", "reminder_surcharge_id"};
    }

    public String getIdentifier() {
        return identifier;
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }
    
    @Override
    public int getId() {
        return id;
    }

    public String getProduct_type() {
        return product_type;
    }

    public int getPayment_type() {
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

    public String getPersonal_template() {
        return personal_template;
    }

    public String getSome_payments() {
        return some_payments;
    }

    public String getShow_in_catalogue() {
        return show_in_catalogue;
    }

    public String getTax() {
        return tax;
    }

    public String getShow_for_partners() {
        return show_for_partners;
    }

    public void setId(int id) {
        this.id = id;
    }

    public void setProduct_type(String product_type) {
        this.product_type = product_type;
    }

    public void setPayment_type(int payment_type) {
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
    
    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public int getGroup_id() {
        return group_id;
    }

    public void setGroup_id(int group_id) {
        this.group_id = group_id;
    }

    public String getPersonal_topic() {
        return personal_topic;
    }

    public void setPersonal_topic(String personal_topic) {
        this.personal_topic = personal_topic;
    }

    public String getPersonal_text() {
        return personal_text;
    }

    public void setPersonal_text(String personal_text) {
        this.personal_text = personal_text;
    }

    public int getReminder_payment_id() {
        return reminder_payment_id;
    }

    public void setReminder_payment_id(int reminder_payment_id) {
        this.reminder_payment_id = reminder_payment_id;
    }

    public int getReminder_surcharge_id() {
        return reminder_surcharge_id;
    }

    public void setReminder_surcharge_id(int reminder_surcharge_id) {
        this.reminder_surcharge_id = reminder_surcharge_id;
    }
    
    
}