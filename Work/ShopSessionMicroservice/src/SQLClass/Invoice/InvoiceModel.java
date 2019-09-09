package SQLClass.Invoice;

import MyInterface.ModelInterface;
import SQLClass.Product.ShortProduct;
import java.sql.Timestamp;
import java.util.ArrayList;

public class InvoiceModel implements ModelInterface {
    int id = -1;
    String email = "";
    String name = ""; 
    String phone = ""; 
    String comment = "";
    String telegram = "";
    String vk = "";
    String facebook = "";
    String postcode = "";
    String country = "";
    String region = "";
    String city = "";
    String address = "";
    String channel = "";
    String source = "";
    String campaign = "";
    String advertising_ads = "";
    String advertising_keys = ""; 
    int user_id = -1;
    String customer_time = ""; 
    String account_validity_period_bool = "";
    double amount = -1.0; 
    Timestamp account_validity_period;
    Timestamp date_create; 
    String payment_status = "";
    String ip_address  = "";
        
    ArrayList<ShortProduct> products = new ArrayList<ShortProduct>();
    
    public InvoiceModel(int id, int user_id, String email, String name, String phone, double amount, String customer_time, Timestamp account_validity_period,
            String account_validity_period_bool, Timestamp date_create, String comment, String telegram, String vk, String facebook, String postcode,
            String country, String region, String city, String address, String channel, String source, String campaign, String advertising_ads, 
            String advertising_keys, String payment_status, String ip_address) {
        this.id = id;
        this.user_id = user_id;
        this.email = email;
        this.name = name;
        this.phone = phone;
        this.amount = amount;
        this.comment = comment;
        this.telegram = telegram;
        this.vk = vk;
        this.facebook = facebook;
        this.postcode = postcode;
        this.country = country;
        this.region = region;
        this.city = city;
        this.address = address;
        this.channel = channel;
        this.source = source;
        this.campaign = campaign;
        this.advertising_ads = advertising_ads;
        this.advertising_keys = advertising_keys;
        this.account_validity_period = account_validity_period; 
        this.account_validity_period_bool = account_validity_period_bool;
        this.date_create = date_create;
        this.customer_time = customer_time;
        this.payment_status = payment_status;
        this.ip_address = ip_address;
    }
    
    public InvoiceModel(int id, int user_id, String email, String name, String phone, double amount, String customer_time, Timestamp account_validity_period,
            String account_validity_period_bool, Timestamp date_create, String comment, String telegram, String vk, String facebook, String postcode,
            String country, String region, String city, String address, String channel, String source, String campaign, String advertising_ads, 
            String advertising_keys, String payment_status, ArrayList<ShortProduct> products, String ip_address) {
        this.id = id;
        this.user_id = user_id;
        this.email = email;
        this.name = name;
        this.phone = phone;
        this.amount = amount;
        this.comment = comment;
        this.telegram = telegram;
        this.vk = vk;
        this.facebook = facebook;
        this.postcode = postcode;
        this.country = country;
        this.region = region;
        this.city = city;
        this.address = address;
        this.channel = channel;
        this.source = source;
        this.campaign = campaign;
        this.advertising_ads = advertising_ads;
        this.advertising_keys = advertising_keys;
        this.account_validity_period = account_validity_period; 
        this.account_validity_period_bool = account_validity_period_bool;
        this.date_create = date_create;
        this.customer_time = customer_time;
        this.payment_status = payment_status;
        this.products = products;
        this.ip_address = ip_address;
    }
        
    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), email, name, phone, String.valueOf(amount), comment, telegram, vk, facebook, postcode, country, region, city,
        address, channel, source, campaign, advertising_ads, advertising_keys, account_validity_period.toString(), account_validity_period_bool, date_create.toString(), 
        customer_time, payment_status, ip_address};
    }
    
    @Override
    public String[] getSqlArray() {
        return new String[] {String.valueOf(id),
            String.valueOf(user_id), 
            "'" + email.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + phone.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            String.valueOf(amount), 
            "'" + comment.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + telegram.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + vk.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + facebook.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + postcode.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + country.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + region.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + city.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + address.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + channel.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + source.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + campaign.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + advertising_ads.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + advertising_keys.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + account_validity_period.toString() + "'", 
            "'" + account_validity_period_bool.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + date_create.toString() + "'", 
            "'" + customer_time.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + payment_status.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"};
    }
    
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "email", "name", "phone", "amount", "comment", "telegram", "vk", "facebook", "postcode", "country", 
            "region", "city", "address", "channel", "source", "campaign", "advertising_ads", "advertising_keys", "account_validity_period", 
            "account_validity_period_bool", "date_create", "customer_time", "payment_status", "ip_address"}; 
    }

    public ArrayList<ShortProduct> getProducts() {
        return products;
    }

    public void setProducts(ArrayList<ShortProduct> products) {
        this.products = products;
    }
    
    public String getPayment_status() {
        return payment_status;
    }

    public void setPayment_status(String payment_status) {
        this.payment_status = payment_status;
    }
    
    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public String getCustomer_time() {
        return customer_time;
    }

    public void setCustomer_time(String customer_time) {
        this.customer_time = customer_time;
    }

    public String getAccount_validity_period_bool() {
        return account_validity_period_bool;
    }

    public void setAccount_validity_period_bool(String account_validity_period_bool) {
        this.account_validity_period_bool = account_validity_period_bool;
    }
    
    @Override
    public int getId() {
        return id;
    }

    public void setId(int Id) {
        this.id = Id;
    }

    public String getEmail() {
        return email;
    }

    public void setEmail(String email) {
        this.email = email;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getPhone() {
        return phone;
    }

    public void setPhone(String phone) {
        this.phone = phone;
    }

    public String getComment() {
        return comment;
    }

    public void setComment(String comment) {
        this.comment = comment;
    }

    public String getTelegram() {
        return telegram;
    }

    public void setTelegram(String telegram) {
        this.telegram = telegram;
    }

    public String getVk() {
        return vk;
    }

    public void setVk(String vk) {
        this.vk = vk;
    }

    public String getFacebook() {
        return facebook;
    }

    public void setFacebook(String facebook) {
        this.facebook = facebook;
    }

    public String getPostcode() {
        return postcode;
    }

    public void setPostcode(String postcode) {
        this.postcode = postcode;
    }

    public String getCountry() {
        return country;
    }

    public void setCountry(String country) {
        this.country = country;
    }

    public String getRegion() {
        return region;
    }

    public void setRegion(String region) {
        this.region = region;
    }

    public String getCity() {
        return city;
    }

    public void setCity(String city) {
        this.city = city;
    }

    public String getAddress() {
        return address;
    }

    public void setAddress(String address) {
        this.address = address;
    }

    public double getAmount() {
        return amount;
    }

    public void setAmount(double amount) {
        this.amount = amount;
    }

    public String getChannel() {
        return channel;
    }

    public void setChannel(String channel) {
        this.channel = channel;
    }

    public String getSource() {
        return source;
    }

    public void setSource(String source) {
        this.source = source;
    }

    public String getCampaign() {
        return campaign;
    }

    public void setCampaign(String campaign) {
        this.campaign = campaign;
    }

    public String getAdvertising_ads() {
        return advertising_ads;
    }

    public void setAdvertising_ads(String advertisingAds) {
        this.advertising_ads = advertisingAds;
    }

    public String getAdvertising_keys() {
        return advertising_keys;
    }

    public void setAdvertisingKeys(String advertising_keys) {
        this.advertising_keys = advertising_keys;
    }
    
    public Timestamp getAccount_validity_period() {
        return account_validity_period;
    }

    public void setAccount_validity_period(Timestamp account_validity_period) {
        this.account_validity_period = account_validity_period;
    }

    public Timestamp getDate_create() {
        return date_create;
    }

    public void setDate_create(Timestamp date_create) {
        this.date_create = date_create;
    }

    public String getIp_address() {
        return ip_address;
    }

    public void setIp_address(String ip_address) {
        this.ip_address = ip_address;
    }
    
}
