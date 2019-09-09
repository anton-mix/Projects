package services;

import java.util.Set;
import javax.ws.rs.core.Application;

@javax.ws.rs.ApplicationPath("webresources")
public class ApplicationConfig extends Application {

    @Override
    public Set<Class<?>> getClasses() {
        Set<Class<?>> resources = new java.util.HashSet<>();
        addRestResourceClasses(resources);
        return resources;
    }

    private void addRestResourceClasses(Set<Class<?>> resources) {
        resources.add(services.ManagerContactsRest.class);
        resources.add(services.ManagerDevLogRest.class);
        resources.add(services.ManagerFilesRest.class);
        resources.add(services.ManagerGroupsContactsRest.class);
        resources.add(services.ManagerGroupsProductsRest.class);
        resources.add(services.ManagerHttpSettingsRest.class);
        resources.add(services.ManagerInvoiceEventsRest.class);
        resources.add(services.ManagerInvoicesRest.class);
        resources.add(services.ManagerMessagesRest.class);
        resources.add(services.ManagerNotificationLogRest.class);
        resources.add(services.ManagerPayments.class);
        resources.add(services.ManagerProductsRest.class);
        resources.add(services.ManagerReminderPaymentsMessageRest.class);
        resources.add(services.ManagerReminderPaymentsRest.class);
        resources.add(services.ManagerReminderSurchargeMessageRest.class);
        resources.add(services.ManagerReminderSurchargeRest.class);
        resources.add(services.ManagerTypesGroupsContactsRest.class);
        resources.add(services.ManagerUsersRest.class);
    }
}
