package Jpetra;
/*
 * Vector.java
 *
 * Created on May 30, 2001, 3:00 PM
 */


/**
 *
 * @author  Mike Heroux
 * @author  Michael William Boldt
 * @version 
 */
public class Vector extends MultiVector {
    
    public Vector(BlockMap map) throws JpetraException {
        super(map, 1);
    }
    
    /** Creates new Vector */
    public Vector(String copyView, BlockMap map, double [][] values)
    throws JpetraException {
        super(copyView, map, values, 1);
    }
    
    public Vector(String copyView, Vector source, int index)
    throws JpetraException {
        super(copyView, source, index, 1);
    }
}
