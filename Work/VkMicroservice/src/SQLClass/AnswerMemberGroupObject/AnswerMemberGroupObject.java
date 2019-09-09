package SQLClass.AnswerMemberGroupObject;

import java.util.ArrayList;


public class AnswerMemberGroupObject {
    int count;
    ArrayList<String> ids;
    public AnswerMemberGroupObject(int count, ArrayList<String> ids){
        this.count = count;
        this.ids = ids;
    }

    public int getCount() {
        return count;
    }

    public void setCount(int count) {
        this.count = count;
    }

    public ArrayList<String> getIds() {
        return ids;
    }

    public void setIds(ArrayList<String> ids) {
        this.ids = ids;
    }
    
    
}
