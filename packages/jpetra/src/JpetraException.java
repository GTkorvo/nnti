/*
 * PetraException.java
 *
 * Created on June 1, 2001, 3:10 PM
 */


package Jpetra;

/**
 *
 * @author  mwboldt
 * @version 
 */
public class JpetraException extends java.lang.Exception {

    /**
     * Creates new <code>PetraException</code> without detail message.
     */
    public JpetraException() {
    }


    /**
     * Constructs an <code>PetraException</code> with the specified detail message.
     * @param msg the detail message.
     */
    public JpetraException(String msg) {
        super(msg);
    }
}


