/*
 * Map.java
 *
 * Created on June 7, 2001, 7:01 PM
 */

package Jpetra;

/**
 *
 * @author  mwboldt
 * @version 
 */
public class Map extends BlockMap {

    /** Creates new Map */
    public Map(int numGlobalElements, int indexBase, Comm comm) 
    //throws JpetraException 
    {
        super(numGlobalElements, 1, indexBase, comm);
    }
    
    public Map(int numGlobalElements, int numProcessElements, int indexBase, Comm comm) 
    //throws JpetraException 
    {
        super(numGlobalElements, numProcessElements, 1, indexBase, comm);
    }
    
    public Map(int numGlobalElements, int numProcessElements, int [] globalElements,
    int indexBase, Comm comm) 
    //throws JpetraException 
    {
        super(numGlobalElements, numProcessElements, globalElements, 1, indexBase, comm);
    }
    
    public Map(Map map) 
    //throws JpetraException 
    {
        super(map);
    }
}
