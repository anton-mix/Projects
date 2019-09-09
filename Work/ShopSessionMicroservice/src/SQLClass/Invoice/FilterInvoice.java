package SQLClass.Invoice;

import MyInterface.ModelInterface;
import SQLClass.Product.ShortProduct;
import java.util.ArrayList;

public class FilterInvoice implements ModelInterface {
    //fild in Invoice BD
    String payment_status = "";
    String channel = "";
    String source = "";
    String campaign = "";
    String advertising_ads = "";
    String advertising_keys = "";
    String postcode = "";
    int user_id = -1;
    
    //option
    int composition = -1;
    
    //other
    String partner = "";
    int click_order = -1;
    String private_manager = "";
    String responsible_manager = "";
    String client_info = "";
    String delivery_status = "";
    String payment_type = "";
    String autopayment = "";   
    String query = "";
    ArrayList<ShortProduct> products;
    int limit = 0;
    int offset = 0;
    
    public String getFullText() {
        return "email, name, phone, payment_status";
    }

    public String[] getArray() {
        return new String[] {payment_status, channel, source, campaign, advertising_ads, advertising_keys, postcode, String.valueOf(user_id),
    String.valueOf(composition), partner, String.valueOf(click_order), private_manager, responsible_manager, client_info, 
    delivery_status, payment_type, autopayment, query};
    }
    
    public String[] getSqlArray() {
        return new String[] {
            "'" + payment_status.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + channel.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + source.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + campaign.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + advertising_ads.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + advertising_keys.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + postcode.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'" ,
            String.valueOf(user_id)};
    }
     
    public String[] getKeyArray() {
        return new String[] {"payment_status", "channel", "source", "campaign", "advertising_ads", "advertising_keys", "postcode", "user_id", "composition",
            "partner" ,"click_order", "private_manager", "responsible_manager", "client_info", "delivery_status", "payment_type", "autopayment", "query"}; 
    }

    public String getQuery() {
        return query;
    }

    public void setQuery(String query) {
        this.query = query;
    }
    
    public String getPayment_status() {
        return payment_status;
    }

    public void setPayment_status(String payment_status) {
        this.payment_status = payment_status;
    }

    public ArrayList<ShortProduct> getProducts() {
        return products;
    }

    public void setProducts(ArrayList<ShortProduct> products) {
        this.products = products;
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

    public void setAdvertising_ads(String advertising_ads) {
        this.advertising_ads = advertising_ads;
    }

    public String getAdvertising_keys() {
        return advertising_keys;
    }

    public void setAdvertising_keys(String advertising_keys) {
        this.advertising_keys = advertising_keys;
    }

    public String getPostcode() {
        return postcode;
    }

    public void setPostcode(String postcode) {
        this.postcode = postcode;
    }

    public int getComposition() {
        return composition;
    }

    public void setComposition(int composition) {
        this.composition = composition;
    }

    public String getPartner() {
        return partner;
    }

    public void setPartner(String partner) {
        this.partner = partner;
    }

    public int getClick_order() {
        return click_order;
    }

    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }
    
    public void setClick_order(int click_order) {
        this.click_order = click_order;
    }

    public String getPrivate_manager() {
        return private_manager;
    }

    public void setPrivate_manager(String private_manager) {
        this.private_manager = private_manager;
    }

    public String getResponsible_manager() {
        return responsible_manager;
    }

    public void setResponsible_manager(String responsible_manager) {
        this.responsible_manager = responsible_manager;
    }

    public String getClient_info() {
        return client_info;
    }

    public void setClient_info(String client_info) {
        this.client_info = client_info;
    }

    public String getDelivery_status() {
        return delivery_status;
    }

    public void setDelivery_status(String delivery_status) {
        this.delivery_status = delivery_status;
    }

    public String getPayment_type() {
        return payment_type;
    }

    public void setPayment_type(String payment_type) {
        this.payment_type = payment_type;
    }

    public String getAutopayment() {
        return autopayment;
    }

    public void setAutopayment(String autopayment) {
        this.autopayment = autopayment;
    }

    public int getLimit() {
        return limit;
    }

    public void setLimit(int limit) {
        this.limit = limit;
    }

    public int getOffset() {
        return offset;
    }

    public void setOffset(int offset) {
        this.offset = offset;
    }
    
    
}
