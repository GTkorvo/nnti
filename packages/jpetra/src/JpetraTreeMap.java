/*
 * JpetraTreeMap.java
 *
 * Created on February 7, 2004, 9:24 PM
 */

package Jpetra;

import java.util.TreeMap;

/**
 *
 * @author  jc
 */
public class JpetraTreeMap extends TreeMap {
    
    /** Creates a new instance of JpetraTreeMap */
    public JpetraTreeMap() {
        super();
    }
    
    public int getInt(int value) {
        return ((Integer) super.get(new Integer(value))).intValue();
    }
}
