/*
 * SerialPlatform.java
 *
 * Created on February 9, 2004, 2:33 PM
 */

package Jpetra;

/**
 *
 * @author  jc
 */
public class SerialPlatform implements Platform {
    
    /** Creates a new instance of SerialPlatform */
    public SerialPlatform() {
    }
    
    public SerialPlatform (SerialPlatform platform) {
        
    }
    
    public Comm createComm() {
        return new SerialComm();
    }
    
    public Directory createDirectory(ElementSpace elementSpace) {
        return new SerialDirectory(elementSpace);
    }
    
    public Distributor createDistributor() {
        return new SerialDistributor();
    }
    
    public Platform clone(Platform platform) {
        return new SerialPlatform((SerialPlatform) platform);
    }
}
