/*
 * VectorSpace.java
 *
 * Created on February 9, 2004, 2:45 PM
 */

package Jpetra;

import java.io.*;

/**
 *
 * @author  jc
 */
public class VectorSpace extends JpetraObject {
    	boolean blockspace;
	int zero;
	int one;
	int indexBase;
	int numMyEntries;
	int numGlobalEntries;
	ElementSpace elementSpace;
	// BlockElementSpace lockElementSpace;
	Platform platform;
	Comm comm;
    
    /** Creates a new instance of VectorSpace */
    public VectorSpace(ElementSpace elementSpace, Platform platform) {
        this.blockspace = false;
        this.zero = 0;
        this.one = 1;
        this.indexBase = elementSpace.getIndexBase();
        this.numMyEntries = elementSpace.getNumMyElements();
        this.numGlobalEntries = elementSpace.getNumGlobalElements();
        
        this.elementSpace = new ElementSpace(elementSpace);
        this.platform = platform.clone(platform);
        this.comm = this.platform.createComm();  // Tpetra had platform.createScalarComm()
    }
    
    /* for BlockElementSpace which isn't implimented yet
    public VectorSpace (BlockElementSpace blockElementSpace, Platform platform) {
        this.blockspace = true;
        this.zero = 0;
        this.one = 1;
        this.indexBase = blockElementSpace.elementSpace().getIndexBase();
        this.numMyEntries = blockElementSpace.getNumMyPoints();
        this.numGlobalEntries = blockElementSpace.getNumGlobalPoints();
        
        this.blockElementSpace = new BlockElementSpace(blockElementSpace);
        this.elementSpace = blockElementSpace.generateCompatibleElementSpace();
        this.platform = new Platform(platform);
        this.comm = this.platform.createComm();
    }*/
    
	
	//! Tpetra::VectorSpace copy constructor.
	public VectorSpace(VectorSpace vectorSpace) {
		this.blockspace = vectorSpace.blockspace;
	this.zero = vectorSpace.zero;
		this.one = vectorSpace.one;
		this.indexBase = vectorSpace.indexBase;
		this.numMyEntries = vectorSpace.numMyEntries;
	this.numGlobalEntries = vectorSpace.numGlobalEntries;
		this.elementSpace = vectorSpace.elementSpace;
		//this.blockElementSpace = vectorSpace.blockElementSpace;
		this.platform = vectorSpace.platform;
		this.comm = vectorSpace.comm;
        }
	
	
	//@{ \name VectorSpace Attribute Methods
	
	//! Returns the number of entries in this VectorSpace.
	public int getNumGlobalEntries() {return(numGlobalEntries);}
	
	//! Returns the number of entries belonging to the calling image.
	public int getNumMyEntries() {return(numMyEntries);}
	
	//! Returns the index base for this VectorSpace.
	public int getIndexBase() {return(indexBase);}
	
	//! Min/Max Indices
	public int getMinLocalIndex() {return(zero);}
	public int getMaxLocalIndex() {return(zero + getNumMyEntries());}
	public int getMinGlobalIndex() {return(getIndexBase());}
	public int getMaxGlobalIndex() {return(getIndexBase() + getNumGlobalEntries());}
	
	//! Return the local index for a given global index
	/*! If this VectorSpace was created using a BlockElementSpace,
		  LIDs and GIDs from the compatible ElementSpace will be used.
	*/
	public int getLocalIndex(int globalIndex) {
			return(elementSpace().getLID(globalIndex));
	}
	
	//! Return the global index for a given local index
	/*! If this VectorSpace was created using a BlockElementSpace,
		  LIDs and GIDs from the compatible ElementSpace will be used.
	*/
	public int getGlobalIndex(int localIndex) {
			return(elementSpace().getGID(localIndex));
	}
	
	//@}
	
	//@{ \name Vector Creation
	
	//! Creates a Tpetra::Vector, with all entries set to 0.
	/*public Vector createVector() {
		Vector vector = new Vector(this);
		return(vector);
	}*/
	
	//@{ \name Boolean Tests
	
	//! Returns true if the VectorSpace passed in is compatible with this VectorSpace.
	public boolean isCompatible(VectorSpace vectorSpace) {
		// first check global length and local lenght on this image
		if(vectorSpace.getNumGlobalEntries() != getNumGlobalEntries() ||
			 vectorSpace.getNumMyEntries() != getNumMyEntries())
			return(false);

		// then check to make sure distribution is the same
		double sameNumLocal = 1.0; // we already know the length matches on this image
		double sameNumGlobal = 0.0;
		if(vectorSpace.getNumMyEntries() == getNumMyEntries())
			sameNumLocal = 1.0;
		double[] temp = {sameNumLocal};
                temp = comm().minAll(temp);
		sameNumGlobal = temp[0];
                return(sameNumGlobal == 1.0);
	}
	
	//! Returns true if the VectorSpace passed in is identical to this VectorSpace. Also implemented through the == and != operators.
	public boolean isSameAs(VectorSpace vectorSpace) {
		/*if(blockspace)
			return(blockElementSpace().isSameAs(vectorSpace.blockElementSpace())); // compare BlockElementSpaces
		else */
			return(elementSpace().isSameAs(vectorSpace.elementSpace())); // compare ElementSpaces
	}
        
	//@}
	
	//@{ \name Misc.
	
	//! Prints the VectorSpace object to the output stream.
	/*! An << operator is inherited from Tpetra::Object, which uses the print method.*/
	public void print(PrintStream os) {
		int myImageID = comm().getVnodeID();
		int numImages = comm().getNumVnodes();
		
		for (int imageCtr = 0; imageCtr < numImages; imageCtr++) {
			if (myImageID == imageCtr) {
				if (myImageID == 0) {
					os.println("\nNumber of Global Entries  = " + getNumGlobalEntries());
					os.println("Index Base                = " + getIndexBase());
				}
				os.println("\nVnodeID = " + myImageID);
				os.println("Number of Local Entries   = " + getNumMyEntries());
				os.println("");
			}
		}
		/*if(blockspace_) {
			os << "Built on a BlockElementSpace" << endl;
			blockElementSpace().print(os);
			os << "Compatible ElementSpace:" << endl;
			elementSpace().print(os);
		}
		else {*/
			os.println("Built on an ElementSpace");
			elementSpace().print(os);
		/*}*/
	};
	
	
	//! Access function for the Tpetra::Platform and Tpetra::Comm communicators.
	public Platform platform() {return(platform);}
	public Comm comm() {return(comm);}
	
	private ElementSpace elementSpace() {return(elementSpace);};
	//private BlockElementSpace blockElementSpace() {return(blockElementSpace);};
	

}
