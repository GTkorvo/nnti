package Jpetra;

import java.io.PrintStream;

/*
 * Class.java
 *
 * Created on September 8, 2003, 2:06 PM
 */

/**
 *
 * @author  Jason Cross
 */
public class Output {
    private String label;
    private boolean rootPrint;
    private boolean vnodesPrint;
    private PrintStream rootStream;
    private PrintStream vnodesStream;
    
    public Output(String label, boolean rootPrint, PrintStream rootStream, boolean vnodesPrint, PrintStream vnodesStream) {
        this.label = label;
        this.rootPrint = rootPrint;
        this.rootStream = rootStream;
        this.vnodesPrint = vnodesPrint;
        this.vnodesStream = vnodesStream;
    }
    
    public void print (String message) {
        if (this.rootPrint) {
            rootStream.print(this.label + message);
        }
        if (this.vnodesPrint) {
            vnodesStream.print(this.label + message);    
        }
    }

    public void println (String message) {
        if (this.rootPrint) {
            rootStream.println(this.label + message);
        }
        if (this.vnodesPrint) {
            vnodesStream.println(this.label + message);    
        }
    }
    
    public boolean isRootPrint () {
        return this.rootPrint;
    }
     
    public boolean isVnodesPrint () {
        return this.vnodesPrint;
    }
    
    public void setRootPrint (boolean rootPrint) {
        this.rootPrint = rootPrint;
    } 
    
    public void setVnodesPrint (boolean vnodesPrint) {
        this.vnodesPrint = vnodesPrint;
    }
    
    public void setRootStream (PrintStream rootStream) {
        this.rootStream = rootStream;
    }
    
    public void setVnodesStream (PrintStream vnodesStream) {
        this.vnodesStream = vnodesStream;
    }
}
