/*
 * ElementSpaceData.java
 *
 * Created on February 7, 2004, 6:38 PM
 */

package Jpetra;

/**
 *
 * @author  jc
 */
public class ElementSpaceData extends JpetraObject {
    Comm comm;
    int numGlobalElements;
    int numMyElements;
    int indexBase;
    int minLID;
    int maxLID;
    int minMyGID;
    int maxMyGID;
    int minAllGID;
    int maxAllGID;
    boolean contiguous;
    boolean global;
    JpetraTreeMap lgMap;
    JpetraTreeMap glMap;
    int[] myGlobalElements;
    Directory directory;
    Platform platform;
    
    /** Creates a new instance of ElementSpaceData */
    public ElementSpaceData(
    int indexBase, int numGlobalElements,int numMyElements,int minAllGID,int maxAllGID,int minMyGID,
    int maxMyGID,
    JpetraTreeMap lgMap,
    JpetraTreeMap glMap,
    boolean contiguous,
    Platform platform,
    Comm comm) {
        this.platform = platform;
        this.comm = comm;
        this.numGlobalElements = numGlobalElements;
        this.numMyElements = numMyElements;
        this.indexBase = indexBase;
        
        // fix		, minLID_(Teuchos::OrdinalTraits<OrdinalType>::zero())
        // fix		, maxLID_(minLID_ + numMyElements_ - Teuchos::OrdinalTraits<OrdinalType>::one())
        this.minMyGID = minMyGID;
        this.maxMyGID = maxMyGID;
        this.minAllGID = minAllGID;
        this.maxAllGID = maxAllGID;
        this.contiguous = contiguous;
        this.global = this.checkGlobalness();
        this.lgMap = lgMap;
        this.glMap = glMap;
        this.myGlobalElements = null;
        this.directory = null;
    }
    
        /*~ElementSpaceData() {
                if(Directory_ != 0) {
                        delete Directory_;
                        Directory_ = 0;
                }
                if(myGlobalElements_ != 0) {
                        delete [] myGlobalElements_;
                        myGlobalElements_ = 0;
                }
                if(Comm_ != 0) {
                        delete Comm_;
                        Comm_ = 0;
                }
        };*/
    
    private boolean checkGlobalness() {
        boolean global = false;
        if(comm.getNumVnodes() > 1) {
            /* huh? assume true 
            int localRep = 0;
            int allLocalRep;
            if(this.numGlobalElements == this.numMyElements)
                localRep = 1;
            comm.minAll(&localRep, &allLocalRep, 1);
            if(allLocalRep != 1)
                global = true;
             */
            global = true;
        }
        return(global);
    }
    
    //! Copy constructor (declared but not defined, do not use)
    //ElementSpaceData(ElementSpaceData<OrdinalType> const& Source);
    //! Assignment operator (declared but not defined, do not use)
    //ElementSpaceData<OrdinalType>& operator = (ElementSpaceData<OrdinalType> const& Source);
    
}
