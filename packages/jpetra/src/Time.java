/*
 * Time.java
 *
 * Created on July 18, 2001, 2:17 PM
 */

package Jpetra;
import java.util.Date;

/**
 *
 * @author  Michael William Boldt
 * @version 
 */
public class Time extends JpetraObject {

    private double startTime;
    private Comm comm;
    private Date date;

    /** Creates new Time */
    public Time(Comm comm) {
        this.comm = comm;
        date = new Date();
        startTime = (double)date.getTime() * 1000.0;
    }

    public Time(Time time) {
        this.comm = time.comm;
        date = new Date();
    }

    void resetStartTime() {
        startTime = (double)date.getTime() * 1000.0;
        return;
    }

    double elapsedTime() {
	    return((double)date.getTime() * 1000.0 - startTime);
    }

}
