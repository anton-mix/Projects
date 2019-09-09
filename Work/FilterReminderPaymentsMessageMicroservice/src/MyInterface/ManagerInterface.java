package MyInterface;

import MyException.MyNotAuthorizedException;
import MyException.MyParseException;
import MyException.MyTypeException;
import SQLClass.SQLResult;
import com.google.gson.Gson;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Calendar;
import java.util.TimeZone;

public interface ManagerInterface {
    
    public default <T extends ModelInterface> String SqlFilter (T item, int Limit, int Offset) throws MyParseException, SQLException, MyTypeException, MyNotAuthorizedException {
        String where_sql = "";
        String[] attribute_invoice = item.getSqlArray();
        String[] key_invoice = item.getKeyArray();
        if("".equals(item.getUser_id())) {
            throw new MyNotAuthorizedException("");
        }
        
        try {
            for (int i = 0; i < attribute_invoice.length; i++) {
                if (!attribute_invoice[i].equals("''") && !attribute_invoice[i].equals("-1") && !attribute_invoice[i].equals("-1.0") && !attribute_invoice[i].equals("'1970-01-02 00:00:00.0'")) {
                    if ("".equals(where_sql)) {
                        where_sql = "WHERE " + key_invoice[i] + "=" + attribute_invoice[i];
                    } else {
                        where_sql = where_sql + " AND " + key_invoice[i] + "=" + attribute_invoice[i] + " ";
                    }
                }
            }
            
            if(item.getQuery().length() != 0 && item.getQuery().replaceAll(" ", "").length() != 0) {
                if(!"".equals(where_sql)) {
                   where_sql = where_sql + " AND MATCH (" + item.getFullText() + ") AGAINST ('" + item.getQuery() + "')"; 
                }
                else {
                    where_sql = "WHERE MATCH (" + item.getFullText() + ") AGAINST ('" + item.getQuery() + "')";         
                }
            }
            
            if (Limit > 0) {
                where_sql = where_sql + " LIMIT " + Limit;
                if (Offset >= 0) {
                    where_sql = where_sql + " OFFSET " + Offset;
                }
            }
        } 
        catch (Exception ex) {
            throw (MyParseException)ex;
        }
         
        return where_sql;
    }
    
    public default String Select (String sql_query, Connection connection_with_db, Gson gson) throws SQLException, MyTypeException {
        Statement statement_with_db = null;
        ResultSet result_with_db = null;
        
        ResultSetMetaData rsmd;
        JsonArray jsonArray = new JsonArray();
        
        try {
            statement_with_db = connection_with_db.createStatement();
            result_with_db = statement_with_db.executeQuery(sql_query);
            rsmd = result_with_db.getMetaData();
            while (result_with_db.next()) {
                JsonObject obj = new JsonObject();
                for (int i = 1; i <= rsmd.getColumnCount(); i++) {
                    switch(rsmd.getColumnType(i)) {
                        case java.sql.Types.INTEGER:
                            obj.addProperty(rsmd.getColumnName(i).toLowerCase(), result_with_db.getInt(i));
                            break;
                        case java.sql.Types.DECIMAL:
                            obj.addProperty(rsmd.getColumnName(i).toLowerCase(), result_with_db.getDouble(i));
                            break;
                        case java.sql.Types.VARCHAR:
                            obj.addProperty(rsmd.getColumnName(i).toLowerCase(), result_with_db.getString(i));
                            break;
                        case java.sql.Types.TIMESTAMP:
                            if(result_with_db.getTimestamp(i, Calendar.getInstance(TimeZone.getDefault())) != null)
                                obj.addProperty(rsmd.getColumnName(i).toLowerCase(), result_with_db.getTimestamp(i,Calendar.getInstance(TimeZone.getDefault())).toString());
                            else 
                                obj.addProperty(rsmd.getColumnName(i).toLowerCase(), "1970-01-02 00:00:00");
                            break;
                        case -5:
                            obj.addProperty(rsmd.getColumnName(i).toLowerCase(), result_with_db.getInt(i));
                            break;
                        default: {
                            throw new MyTypeException("Unknown type " + rsmd.getColumnType(i));
                        }
                    }
                }
                jsonArray.add(obj);
            }

            result_with_db.close();
            statement_with_db.close();
             
            return gson.toJson(jsonArray);

        } 
        catch (SQLException ex) {
            throw ex;
        }
    }
    
    public default <T extends ModelInterface> String ExecuteQuery(String query, Connection connection_with_db, Gson gson) throws SQLException, MyParseException, MyNotAuthorizedException {
        Statement StatementWithDB = null;
        ResultSet result_with_db = null;
               
        try {
            StatementWithDB = connection_with_db.createStatement();
            StatementWithDB.executeUpdate(query);
            StatementWithDB.close();
            
            return gson.toJson(new SQLResult("Ok"));
        }
        catch(SQLException ex)  {
            throw ex;
        }
    }
    
    public default <T extends ModelInterface> String Add(T item, String table, Connection connection_with_db, Gson gson, boolean use_user_id) throws SQLException, MyParseException, MyNotAuthorizedException {
        Statement StatementWithDB = null;
        ResultSet result_with_db = null;
        
        String[] AttributeInvoice = item.getSqlArray();
        String[] KeyInvoice = item.getKeyArray();
        String NeedAttributeInvoice = "", NeedKeyInvoice = "";
        
        if(use_user_id && "".equals(item.getUser_id())) {
            throw new MyNotAuthorizedException("");
        }
        
        try {
            for (int i = 0; i < AttributeInvoice.length; i++) {
                if (!AttributeInvoice[i].equals("''") && !AttributeInvoice[i].equals("-1") && !AttributeInvoice[i].equals("-1.0") && !AttributeInvoice[i].equals("'1970-01-02 00:00:00.0'")) {
                    if ("".equals(NeedAttributeInvoice)) {
                        NeedAttributeInvoice = AttributeInvoice[i];
                        NeedKeyInvoice = KeyInvoice[i];
                    } else {
                        NeedAttributeInvoice = NeedAttributeInvoice + "," +  AttributeInvoice[i];
                        NeedKeyInvoice = NeedKeyInvoice + "," + KeyInvoice[i];
                    }
                }
            }
        } 
        catch (Exception ex) {
            throw (MyParseException)ex;
        }
              
        try {
            StatementWithDB = connection_with_db.createStatement();
            StatementWithDB.executeUpdate("INSERT INTO " + table + " ( " + NeedKeyInvoice + " ) VALUES (" + NeedAttributeInvoice + ");", Statement.RETURN_GENERATED_KEYS);
            result_with_db = StatementWithDB.getGeneratedKeys();
            result_with_db.next();
            int buff = result_with_db.getInt(1);
            
          
            StatementWithDB.close();
            
            return gson.toJson(new SQLResult("Ok",buff, item.getUser_id()));
        }
        catch(SQLException ex)  {
            throw ex;
        }
    }
    
    public default <T extends ModelInterface> String Delete(T item, String table, Connection connection_with_db, Gson gson) throws SQLException, MyNotAuthorizedException {
        Statement StatementWithDB = null;
        
        if("".equals(item.getUser_id())) {
            throw new MyNotAuthorizedException("");
        }
        
        try {
            StatementWithDB = connection_with_db.createStatement();
            StatementWithDB.executeUpdate("DELETE FROM " + table + " WHERE id=" + item.getId() + " AND user_id = '" + item.getUser_id() + "';");
            StatementWithDB.close();
            
            return gson.toJson(new SQLResult("Ok"));
        }
        catch(SQLException ex)  {
            throw ex;
        }
    }
    
    public default <T extends ModelInterface> String Update(T item, String table, Connection connection_with_db, Gson gson) throws SQLException, MyParseException, MyNotAuthorizedException {
        Statement statement_with_db = null;
        String[] attribute_invoice = item.getSqlArray();
        String[] key_invoice = item.getKeyArray();
        String set_string = "" ;
        
        if("".equals(item.getUser_id())) {
            throw new MyNotAuthorizedException("");
        }
        
        try {
            for (int i = 1; i < attribute_invoice.length; i++) {
                if (!attribute_invoice[i].equals("''") && !attribute_invoice[i].equals("-1") && !attribute_invoice[i].equals("-1.0") && !attribute_invoice[i].equals("'1970-01-02 00:00:00.0'")) {
                    if ("".equals(set_string)) {
                        set_string = set_string + key_invoice[i] + "=" + attribute_invoice[i];
                    } else {
                        set_string = set_string + ", " + key_invoice[i] + "=" + attribute_invoice[i] + " ";
                    }
                }
            }
        } 
        catch (Exception ex) {
            throw (MyParseException)ex;
        }
        
        try {
            statement_with_db = connection_with_db.createStatement();
            statement_with_db.executeUpdate("UPDATE " + table + " SET " + set_string + " WHERE id=" + item.getId() + " AND user_id = '" + item.getUser_id() + "';");
            statement_with_db.close();            
            return gson.toJson(new SQLResult("Ok"));
        }
        catch(SQLException ex)  {
            throw ex;
        }
    }
}
