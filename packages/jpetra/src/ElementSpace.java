/*
 * ElementSpace.java
 *
 * Created on February 7, 2004, 6:30 PM
 */

package Jpetra;

import java.util.*;
import java.io.*;
/**
 *
 * @author  jc
 */
public class ElementSpace extends JpetraObject {
    ElementSpaceData elementSpaceData;
    int zero = 0;
    Platform platform;
    
    /** Creates a new instance of ElementSpace */
    public ElementSpace(int numGlobalElements, int indexBase, Platform platform) {
        //const OrdinalType one = Teuchos::OrdinalTraits<OrdinalType>::one();
        //const OrdinalType zero = Teuchos::OrdinalTraits<OrdinalType>::zero();
        int zero = 0;
        int one = 1;
        // initial throws
        if (numGlobalElements < zero) {}
        //throw reportError("numGlobalElements = " + toString(numGlobalElements) + ".  Should be >= " + toString(zero) + ".", -1);
        
        // platform & comm setup
        this.platform = platform;
        Comm comm = platform.createComm();
        int numImages = comm.getNumVnodes();
        int myImageID = comm.getVnodeID();
        
        // compute numMyElements
        int numMyElements = numGlobalElements / numImages;
        int remainder = numGlobalElements % numImages;
        int start_index = myImageID * (numMyElements + one);
        if (myImageID < remainder)
            numMyElements++;
        else
            start_index -= (myImageID - remainder);
        
        // setup lgmap & glmap
        JpetraTreeMap lgMap = new JpetraTreeMap();
        JpetraTreeMap glMap = new JpetraTreeMap();
        
        // setup min/maxs
        int minAllGID = indexBase;
        int maxAllGID = minAllGID + numGlobalElements - one;
        int minMyGID = start_index + indexBase;
        int maxMyGID = minMyGID + numMyElements - one;
        
        // call ESData constructor
        elementSpaceData = new ElementSpaceData(indexBase, numGlobalElements, numMyElements, minAllGID, maxAllGID,
        minMyGID, maxMyGID, lgMap, glMap, true, platform, comm);
        
        // initialize directory
        directorySetup();
    }
    
    public ElementSpace(ElementSpace elementSpace) {
        this.elementSpaceData = elementSpace.elementSpaceData;
        this.platform = elementSpace.platform;
    }
    
    Comm getComm() {return elementSpaceData.comm;}
    /* Local/Global ID Accessor Methods */
    //! Returns the image IDs and corresponding local IDs for a given list of global IDs.
    void getRemoteIDList(int numIDs, int[] GIDList, int[] imageIDList, int[] LIDList)
    {elementSpaceData.directory.getDirectoryEntries(numIDs, GIDList, imageIDList, LIDList);};
    
    //! Returns local ID of global ID passed in, throws exception -1 if not found on this image.
    int getLID(int GID) {
        if(!isMyGID(GID))
            this.println("ERR", "Global ID " + GID + " was not found on this image.");
        else if(isContiguous())
            return(GID - getMinMyGID()); //compute with offset
        else {
            return(elementSpaceData.glMap.getInt(GID));
        }
        
        this.println("ERR", "Shouldn't have hit this point in ElementSpace.getLID().");
        return 0;
    }
    
    //! Returns global ID of local ID passed in, throws exception -1 if not found on this image.
    int getGID(int LID) {
        if(!isMyLID(LID))
            this.println("ERR", "Local ID " + LID + " was not found on this image. (ElementSpace.getLID())");
        else if(isContiguous())
            return(LID + getMinMyGID()); //compute with offset
        else {
            return(elementSpaceData.lgMap.getInt(LID));
        }
        
        this.println("ERR", "Shouldn't have hit this point in ElementSpace.getGID().");
        return 0;
    }
    
    //! Returns true if global ID passed in belongs to the calling image, returns false if it doesn't.
    boolean isMyGID(int GID) {
        if(GID < getMinMyGID() || GID > getMaxMyGID())
            return(false);
        
        if(isContiguous()) return(true);
        
        return (elementSpaceData.glMap.get(new Integer(GID)) != null);
    }
    
    //! Returns true if the local ID passed in belongs to the calling image, returns false if it doesn't.
    boolean isMyLID(int LID) {
        if(LID < getMinLID() || LID > getMaxLID())
            return(false);
        else if(isContiguous())
            return(true);
        else {
            return (elementSpaceData.lgMap.get(new Integer(LID)) != null);
        }
    }
    
    //! Returns the minimum global ID in this ElementSpace.
    int getMinAllGID() {return(elementSpaceData.minAllGID);}
    
    //! Returns the maximum global ID in this ElementSpace.
    int getMaxAllGID() {return(elementSpaceData.maxAllGID);}
    
    //! Returns the minimum global ID owned by this image.
    int getMinMyGID() {return(elementSpaceData.minMyGID);}
    
    //! Returns the maximum global ID owned by this image.
    int getMaxMyGID() {return(elementSpaceData.maxMyGID);}
    
    //! Returns the minimum local ID on the calling image.
    int getMinLID() {return(elementSpaceData.minLID);}
    
    //! Returns the maximum local ID on the calling image.
    int getMaxLID() {return(elementSpaceData.maxLID);}
    
    
    /* Size & Dimension Accessor Methods */
    
    //! Returns the number of elements in this ElementSpace.
    int getNumGlobalElements() {return(elementSpaceData.numGlobalElements);}
    
    //! Returns the number of elements belonging to the calling image.
    int getNumMyElements() {return(elementSpaceData.numMyElements);}
    
    //! Puts list of global elements on this image into the user-provided array.
    void getMyGlobalElements(int[] elementList) {
        if(elementList == null)
            //throw reportError("Pointer does not have child allocated.", 3);
            this.println("ERR", "elementList is null and shouldn't be. (ElementSpace.getMyGlobalElements)");
        else if(isContiguous()) {
            int nME = getNumMyElements();
            int minMyGID = getMinMyGID();
            for(int i = zero; i < nME; i++)
                elementList[i] = minMyGID + i;
        }
        else { // not contiguous
            Iterator iterator = elementSpaceData.lgMap.values().iterator();
            int i=0;
            while(iterator.hasNext()) {
                elementList[i++] = ((Integer) iterator.next()).intValue();
            }
        }
    }
    
    int[] getMyGlobalElements() {
        int nME = getNumMyElements();
        if((elementSpaceData.myGlobalElements == null) && (nME > zero)) {
            elementSpaceData.myGlobalElements = new int[nME];
            getMyGlobalElements(elementSpaceData.myGlobalElements);
        }
        return(elementSpaceData.myGlobalElements);
    }
    
    //! Returns the Index base for this ElementSpace. Normally 0 for C/C++ or 1 for Fortran, but can be anything.
    int getIndexBase() {return(elementSpaceData.indexBase);}
    
    //@}
    
    /*some stuff*/
    private void directorySetup() {
        if(getNumGlobalElements() != this.zero)
            if(elementSpaceData.directory == null)
                elementSpaceData.directory = platform.createDirectory(this); // Make directory
    }
    
    
    /* \name Misc. Boolean Tests */
    
    //! Returns true if this ElementSpace is distributed contiguously, returns false otherwise.
    boolean isContiguous() {return(elementSpaceData.contiguous);};
    
    //! Returns true if this ElementSpace is distributed across more than one image, returns false otherwise.
    boolean isGlobal() {return(elementSpaceData.global);};
    
    //! Returns true if the ElementSpace passed in is identical to this ElementSpace. Also implemented through the == and != operators.
    boolean isSameAs(ElementSpace elementSpace) {
        // Quickest test: See if we share an inner data class
        //if(ElementSpaceData_.shares_resource(ElementSpace.ElementSpaceData_))
        //	return(true);
        
        // Next check other global properties that are easy global attributes
        if(getMinAllGID() != elementSpace.getMinAllGID() ||
        getMaxAllGID() != elementSpace.getMaxAllGID() ||
        getNumGlobalElements() != elementSpace.getNumGlobalElements() ||
        getIndexBase() != elementSpace.getIndexBase() ||
        isGlobal() != elementSpace.isGlobal() ||
        isContiguous() != elementSpace.isContiguous())
            return(false);
        
        // If we get this far, we need to check local properties and then check across
        // all images to see if local properties are all true
        
        // check that maps have the same number of local elements
        int mySameSpace = 1;
        if(getNumMyElements() != elementSpace.getNumMyElements())
            mySameSpace=0;
        
        // then check that GIDs are the same, by checking that the maps match
        // *possible optimization* is this only necessary in a non-contiguous elementspace?
        if(mySameSpace == 1)
            if(elementSpaceData.lgMap != elementSpace.elementSpaceData.lgMap)
                mySameSpace=0;
        
        // Now get min of mySameSpace across all images
        int globalSameSpace = 0;
        int[] tmp = {mySameSpace};
        tmp = getComm().minAll(tmp);
        globalSameSpace = tmp[0];
        // if globalSameSpace is 1, that means none of the images set it to 0,
        // so the ElementSpaces are identical on all images. If this is the case, then we should return true.
        return(globalSameSpace == 1);
    }
    
    //@}
    
    void print(PrintStream os) {
        int zero = 0;
        
        int[] myGlobalElements1 = getMyGlobalElements();
        int myImageID = getComm().getVnodeID();
        int numImages = getComm().getNumVnodes();
        int minLID = getMinLID();
        int nME = getNumMyElements();
        
        for (int imageCtr = zero; imageCtr < numImages; imageCtr++) {
            if (myImageID == imageCtr) {
                if (myImageID == zero) {
                    os.println("\nNumber of Global Elements  = " + getNumGlobalElements());
                    os.println("Maximum of all GIDs        = " + getMaxAllGID());
                    os.println("Minimum of all GIDs        = " + getMinAllGID());
                    os.println("Index Base                 = " + getIndexBase());
                }
                os.println("");
                
                os.println("Number of Local Elements   = " + getNumMyElements());
                os.println("Maximum of my GIDs         = " + getMaxMyGID());
                os.println("Minimum of my GIDs         = " + getMinMyGID());
                os.println("");
                
                os.print("   ImageID" + "    ");
                os.print("       Local Index " + " ");
                os.print("      Global Index " + " ");
                os.println("");
                
                for (int i = zero, lid = minLID; i < nME; i++, lid++) {
                    os.print(myImageID + "    ");
                    os.print(lid + "    ");
                    os.print(myGlobalElements1[i] + "    ");
                    os.println("");
                }
                
            }
            // Do a few global ops to give I/O a chance to complete
            getComm().barrier();
            getComm().barrier();
            getComm().barrier();
        }
    }
    
    
}
