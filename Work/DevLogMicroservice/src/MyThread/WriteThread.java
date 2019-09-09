package MyThread;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class WriteThread implements Runnable {
    private String ReportsAll = "";
    private String ReportsString = "";
    File file;

    public WriteThread() {
        System.out.println("Start MyThread");
    }
    
    @Override
    public void run() {
        while (true) {
            try {
                file = new File("reportsDir");
                if(!file.exists()) {
                    file.mkdir();
                }
                
                System.out.println("Writing....");
                if(!"".equals(ReportsAll)) {
                    doWork(file.getPath(), ReportsAll);
                    ReportsAll = "";
                }
                Thread.sleep(30000);
            } catch (InterruptedException ex) {
                System.out.println("Thread sleep " + ex.toString());
            }
        }
    }
    
    private static void doWork(String Path, String message) {
        FileWriter writer = null;
        try {
            writer = new FileWriter(Path + "/reports.txt", true);
            writer.write(message);
            writer.flush();
            writer.close();
        } catch (IOException ex) {
            System.out.println("error" + ex.getMessage());
        }
    }
    
    public String getReportsString() {
        return ReportsString;
    }

    public void setReportsString(String ReportsString) {
        this.ReportsString = ReportsString;
        ReportsAll = ReportsAll + "\n" + ReportsString;
    }
}
