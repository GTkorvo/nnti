package test;

import Jpetra.*;

/**
 *
 * @author  Jason Cross
 */
public interface JpetraServerApplication {

    public void setComm(Comm comm);
    
    public void runApplication(String[] args);
}
