package MyException;

public class MyNotAuthorizedException extends Exception {
    public MyNotAuthorizedException(String message) {
        super(message);
    }
}
