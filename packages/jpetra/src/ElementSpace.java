package Jpetra;

/**
 *
 * @author  jacross
 */
public class ElementSpace {
    private int numGlobalElements;
    private int[] myGlobalElements;
    private Comm comm;
    
    /** 
      * Constructs an <code>ElementSpace</code> automatically based on the even
      * distribution of global elements to each vnode.
      */
    public ElementSpace(int numGlobalElements, Comm comm) {
        // for both Serial and Parallel
        this.comm = comm;
        this.numGlobalElements = numGlobalElements;
        
        // for serial only
        if (comm.getIsSerial()) {
            this.myGlobalElements = new int[numGlobalElements];
            for (int i=0; i < numGlobalElements; i++) {
                this.myGlobalElements[i] = i;
            }
            
            return;
        }
        
        // for parallel only
        int numIndicesPerVnode = numGlobalElements / comm.getNumVnodes();
        int numRemainderIndices = numGlobalElements % comm.getNumVnodes();
        int numMyGlobalIds = numIndicesPerVnode;
        int numRemainderIndicesToAdd = 0;
        // if the ElementSpace does not map evenly onto all vectors
        // then give vnodes < remainder an additional element
        if (comm.getVnodeID() < numRemainderIndices) {
            numMyGlobalIds++;
            // accounts for the indices owned by vnodes < myVnode
            numRemainderIndicesToAdd = comm.getVnodeID();
        }
        else {
            numRemainderIndicesToAdd = numRemainderIndices;
        }
        myGlobalElements = new int[numMyGlobalIds];
        int myGlobalIdIndex = comm.getVnodeID() * numIndicesPerVnode + numRemainderIndices;
        for (int i=0; i < numMyGlobalIds; i++) {
            myGlobalElements[i] = myGlobalIdIndex++;
        }
        
    }
    
    /**
      * Constructs a user defined <code>ElementSpace</code>.
      */
    public ElementSpace(int[] myGlobalElements, Comm comm) {
        // for both serial and parallel
        this.myGlobalElements = myGlobalElements;
        this.comm = comm;
        
        // serial only
        if (comm.getIsSerial()) {
            this.numGlobalElements = myGlobalElements.length;
            return;
        }
        
        // parallel only
        // need to get the number of global elements by doing a gather of the number
        // of elements on each vnode
        // this.numGlobalElements = ?
    }
    
    public int getNumGlobalElements() {
        return numGlobalElements;
    }
    
    public int getNumMyElements() {
        return myGlobalElements.length;
    }
}
