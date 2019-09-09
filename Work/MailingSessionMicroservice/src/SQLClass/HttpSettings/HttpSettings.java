
package SQLClass.HttpSettings;

import MyInterface.ModelInterface;
import java.util.ArrayList;


public class HttpSettings implements ModelInterface {
    int id = -1; 
    int user_id = -1; 
    
    String emails = "";
    String url_sub = "";
    String use_default_act = "";
    String url_act = "";
    String use_default_whitelist = "";
    String url_unsub = "";
    String use_default_panel = "";
    String auto_unsub = "";
    
    String logo = ""; 
    String logo_link = ""; 
    String description = ""; 
    String contact_icon = ""; 
    String main_contact = ""; 
    String main_contact_link = ""; 
    String additional_contact = ""; 
    String additional_contact_link = ""; 
    String picture = ""; 
    String greeting = ""; 
    String template = ""; 
    String social_message = ""; 
    ArrayList<SocialIcons> social_icons;
    ArrayList<UserEmails> user_emails;
    
    public HttpSettings(int user_id){
        this.id = -1;
        this.user_id = user_id;
        
        this.emails = "";
        this.url_sub = "";
        this.use_default_act = "";
        this.url_act = "";
        this.use_default_whitelist = "";
        this.url_unsub = "";
        this.use_default_panel = "";
        this.auto_unsub = "";
        
        this.logo = "";
        this.logo_link = "";
        this.description = "";
        this.contact_icon = "";
        this.main_contact = "";
        this.main_contact_link = "";
        this.additional_contact = "";
        this.additional_contact_link = "";
        this.picture = "";
        this.greeting = "";
        this.template = "";
        this.social_message = "";
        
        this.social_icons = new ArrayList<SocialIcons>();
        this.user_emails = new ArrayList<UserEmails>();
    }
    
    
    public HttpSettings(int id, int user_id, String emails, String url_sub, String use_default_act,
        String url_act, String use_default_whitelist, String url_unsub, String use_default_panel, String auto_unsub,
        String logo, String logo_link, String description, String contact_icon, String main_contact, 
        String main_contact_link, String additional_contact, String additional_contact_link,
        String picture, String greeting, String template, String social_message, ArrayList<SocialIcons> social_icons,
        ArrayList<UserEmails> user_emails){
        this.id = id;
        this.user_id = user_id;
        
        this.emails = emails;
        this.url_sub = url_sub;
        this.use_default_act = use_default_act;
        this.url_act = url_act;
        this.use_default_whitelist = use_default_whitelist;
        this.url_unsub = url_unsub;
        this.use_default_panel = use_default_panel;
        this.auto_unsub = auto_unsub;
        
        this.logo = logo;
        this.logo_link = logo_link;
        this.description = description;
        this.contact_icon = contact_icon;
        this.main_contact = main_contact;
        this.main_contact_link = main_contact_link;
        this.additional_contact = additional_contact;
        this.additional_contact_link = additional_contact_link;
        this.picture = picture;
        this.greeting = greeting;
        this.template = template;
        this.social_message = social_message;
        
        this.social_icons = social_icons;
        this.user_emails = user_emails;
    }

    @Override
    public String[] getArray() {
        return new String[] {String.valueOf(id), String.valueOf(user_id), emails, url_sub, use_default_act,
            url_act, use_default_whitelist, url_unsub, use_default_panel, auto_unsub, logo, logo_link, description, 
            contact_icon, main_contact, main_contact_link, additional_contact, additional_contact_link, picture,
            greeting, template, social_message};
    }

    @Override
    public String[] getSqlArray() {
        return new String[] {
            String.valueOf(id), 
            String.valueOf(user_id),
            
            "'" + emails.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + url_sub.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + use_default_act.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + url_act.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + use_default_whitelist.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + url_unsub.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + use_default_panel.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + auto_unsub.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            
            "'" + logo.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + logo_link.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + description.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + contact_icon.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + main_contact.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + main_contact_link.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + additional_contact.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + additional_contact_link.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + picture.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + greeting.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'", 
            "'" + template.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + social_message.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
        };
    }
    
    @Override
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "emails", "url_sub", "use_default_act",
            "url_act", "use_default_whitelist", "url_unsub", "use_default_panel", "auto_unsub", 
            "logo", "logo_link", "description", "contact_icon", "main_contact", "main_contact_link", 
            "additional_contact", "additional_contact_link", "picture", "greeting", "template", "social_message"}; 
    }
    
    @Override
    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    @Override
    public int getUser_id() {
        return user_id;
    }

    public void setUser_id(int user_id) {
        this.user_id = user_id;
    }

    public String getLogo() {
        return logo;
    }

    public void setLogo(String logo) {
        this.logo = logo;
    }

    public String getLogo_link() {
        return logo_link;
    }

    public void setLogo_link(String logo_link) {
        this.logo_link = logo_link;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getMain_contact() {
        return main_contact;
    }

    public void setMain_contact(String main_contact) {
        this.main_contact = main_contact;
    }

    public String getMain_contact_link() {
        return main_contact_link;
    }

    public void setMain_contact_link(String main_contact_link) {
        this.main_contact_link = main_contact_link;
    }

    public String getAdditional_contact() {
        return additional_contact;
    }

    public void setAdditional_contact(String additional_contact) {
        this.additional_contact = additional_contact;
    }

    public String getAdditional_contact_link() {
        return additional_contact_link;
    }

    public void setAdditional_contact_link(String additional_contact_link) {
        this.additional_contact_link = additional_contact_link;
    }

    public String getGreeting() {
        return greeting;
    }

    public void setGreeting(String greeting) {
        this.greeting = greeting;
    }

    public String getTemplate() {
        return template;
    }

    public void setTemplate(String template) {
        this.template = template;
    }

    public String getContact_icon() {
        return contact_icon;
    }

    public void setContact_icon(String contact_icon) {
        this.contact_icon = contact_icon;
    }

    public String getPicture() {
        return picture;
    }

    public void setPicture(String picture) {
        this.picture = picture;
    }

    public String getSocial_message() {
        return social_message;
    }

    public void setSocial_message(String social_message) {
        this.social_message = social_message;
    }

    public ArrayList<SocialIcons> getSocial_icons() {
        return social_icons;
    }

    public void setSocial_icons(ArrayList<SocialIcons> social_icons) {
        this.social_icons = social_icons;
    }

    public String getEmails() {
        return emails;
    }

    public void setEmails(String emails) {
        this.emails = emails;
    }

    public String getUrl_sub() {
        return url_sub;
    }

    public void setUrl_sub(String url_sub) {
        this.url_sub = url_sub;
    }

    public String getUse_default_act() {
        return use_default_act;
    }

    public void setUse_default_act(String use_default_act) {
        this.use_default_act = use_default_act;
    }

    public String getUrl_act() {
        return url_act;
    }

    public void setUrl_act(String url_act) {
        this.url_act = url_act;
    }

    public String getUse_default_whitelist() {
        return use_default_whitelist;
    }

    public void setUse_default_whitelist(String use_default_whitelist) {
        this.use_default_whitelist = use_default_whitelist;
    }

    public String getUrl_unsub() {
        return url_unsub;
    }

    public void setUrl_unsub(String url_unsub) {
        this.url_unsub = url_unsub;
    }

    public String getUse_default_panel() {
        return use_default_panel;
    }

    public void setUse_default_panel(String use_default_panel) {
        this.use_default_panel = use_default_panel;
    }

    public String getAuto_unsub() {
        return auto_unsub;
    }

    public void setAuto_unsub(String auto_unsub) {
        this.auto_unsub = auto_unsub;
    }

    public ArrayList<UserEmails> getUser_emails() {
        return user_emails;
    }

    public void setUser_emails(ArrayList<UserEmails> user_emails) {
        this.user_emails = user_emails;
    }

    

    
    
    
    
}
