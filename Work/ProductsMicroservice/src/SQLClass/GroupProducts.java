package SQLClass;

import java.util.List;

public class GroupProducts {
    transient int id = -1;
    transient String user_id = "";
    transient String group_name = "";

    enum GroupProducts_enum {id, user_id, group_name}
    
    public GroupProducts (int id, String user_id, String group_name)  {
        this.id = id;
        this.user_id = user_id;
        this.group_name = group_name;
    }
    
    public GroupProducts (String user_id, String group_name)  {
        this.user_id = user_id;
        this.group_name = group_name;
    }
    
    public GroupProducts (String json) {
        List<String[]> jsonString = GlobalVariables.Parse_json(json);

        for(int i = 0; i < jsonString.size(); i ++) {

            GroupProducts_enum groupProducts_enum = GroupProducts_enum.valueOf(jsonString.get(i)[0]);
            switch(groupProducts_enum) {
                case id:
                    this.id = Integer.valueOf(jsonString.get(i)[1]);
                    break;
                case user_id:
                    this.user_id = jsonString.get(i)[1];
                    break;
                case group_name:
                    this.group_name = jsonString.get(i)[1];
                    break;
            }
        }
    }
    
    public String[] getArray() {
        return new String[] {String.valueOf(id), user_id, group_name};
    }
    
    public String[] getSqlArray() {
        return new String[] {String.valueOf(id),
            "'" + user_id.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'",
            "'" + group_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'"  };
    }
    
    public String getSqlStringArray() {
        return "null ," +
            "'" + user_id.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "'," +
            "'" + group_name.replaceAll("'", "\\\\'").replaceAll("\"", "\\\\\"") + "',";
    }
    
    public String[] getKeyArray() {
        return new String[] {"id", "user_id", "group_name"};
    }
    
    public String getKeyString() {
        return "id, user_id, group_name"; 
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

    public void setId(int id) {
        this.id = id;
    }

    public String getUser_id() {
        return user_id;
    }

    public void setUser_id(String user_id) {
        this.user_id = user_id;
    }

    public String getGroup_name() {
        return group_name;
    }

    public void setGroup_name(String group_name) {
        this.group_name = group_name;
    }

}
