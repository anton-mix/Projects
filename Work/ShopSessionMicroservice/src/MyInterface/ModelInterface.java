package MyInterface;

public interface ModelInterface {
    public String[] getArray();
    public String[] getSqlArray();
    public String[] getKeyArray();
    
    public default int getId() {
        return -1;
    }
    
    public default int getUser_id()  {
        return -1;
    }
    
    public default String getQuery() {
        return "";
    }
    
    public default String getFullText() {
        return "";
    }
    
    public default int getLimit() {
        return 0;
    }
    
    public default int getOffset() {
        return 0;
    }
}
