package Jpetra;
/*
 * BlockMap.java
 *
 * Created on May 29, 2001, 2:21 PM
 */

/**
 * @author  Mike Heroux
 * @author  Michael William Boldt
 * @version 
 */
public class BlockMap extends JpetraObject {
    
    protected int numVnodeElements = 0;
    private int numGlobalElements = 0;
    private int [] globalElements = null;
    private int [] firstElementEntryList = null;
    private int [] equationToBlockList = null;
    private int numVnodeEquations = 0;
    private int numGlobalEquations = 0;
    private int [] processLID = null;
    private int minAllGID = 0;
    private int maxAllGID = 0;
    private int minVnodeGID = 0;
    private int maxVnodeGID = 0;
    private int minLID = 0;
    private int maxLID = 0;
    private int maxVnodeElementSize = 0;
    private int minVnodeElementSize = 0;
    private int maxElementSize = 0;
    private int minElementSize = 0;
    private int indexBase = 0;
    private int elementSize = 0;
    private int [] elementSizeList = null;
    protected Comm comm;
    private Directory directory = null;
    private boolean hasConstantElementSize = false;
    private boolean isLinearMap = false;
    protected boolean isDistributedGlobal = false;
        
    /**
     * Constructor for a Petra-defined uniform linear distribution of
     * constant block size elements.
     *
     * @param numGlobalElements the number of elements to distribute
     * @param elementSize       the number of equations or vectors per element
     * @param indexBase         the minimum index value used for arrays that
     *                          use this map, typically 0 for C/C++ and 1 for Fortran
     * @param comm              the communicator containing information on the number
     *                          of processes
     */
    public BlockMap(int numGlobalElements, int elementSize, int indexBase, Comm comm)
	//throws JpetraException
    {
        this.numGlobalElements = numGlobalElements;
        this.elementSize = elementSize;
        this.indexBase = indexBase;
        this.comm = comm;
        minVnodeElementSize = elementSize;
        maxVnodeElementSize = elementSize;
        minElementSize = elementSize;
        maxElementSize = elementSize;
        hasConstantElementSize = true;
        isLinearMap = true;
        isDistributedGlobal = true;
        
        // if (this.numGlobalElements < 0) throw new JpetraException("numGlobalElements < 0.");
        // if (this.elementSize <= 0) throw new JpetraException("elementSize <= 0.");
	if(this.numGlobalElements < 0) {
	    System.out.println("BlockMap constructor: numGlobalElements < 0");
	    System.exit(1);
	}
	if(this.elementSize <= 0) {
	    System.out.println("BlockMap constructor: elementSize <= 0");
	    System.exit(1);
	}
        
        int numProc = comm.getNumVnodes();
        int myPID = comm.getVnodeID();
        numVnodeElements = this.numGlobalElements / numProc;
        int remainder = this.numGlobalElements % numProc;
        int startIndex = myPID * (numVnodeElements + 1);
        
        if (myPID < remainder) numVnodeElements++;
        else startIndex -= (myPID - remainder);
        
        numGlobalEquations = this.numGlobalElements * this.elementSize;
        numVnodeEquations = numVnodeElements * this.elementSize;
        
        minVnodeElementSize = this.elementSize;
        maxVnodeElementSize = this.elementSize;
        minElementSize = this.elementSize;
        maxElementSize = this.elementSize;
        
        minAllGID = this.indexBase;
        maxAllGID = minAllGID + this.numGlobalElements - 1;
        minVnodeGID = startIndex + this.indexBase;
        maxVnodeGID = minVnodeGID + this.numVnodeElements - 1;
        minLID = 0;
        maxLID = minLID + this.numVnodeElements - 1;
        if (numProc == 1) isDistributedGlobal = false;
        globalToLocalSetup();
    }
    
    /**
     * Constructor for a user-defined linear distribution of constant block size elements.
     *
     * @param numGlobalElements the number of elements to distribute
     * @param numVnodeElements     the number of elements owned by the calling process
     * @param elementSize       the number of equations or vector entries per element
     * @param indexBase         the minimum index value used for arrays that use this
     *                          map, typically 0 for C/C++ and 1 for Fortran
     * @param comm              the communicator containing information on the number
     *                          of processes
     */
    public BlockMap(int numGlobalElements, int numVnodeElements, int elementSize, 
    int indexBase, Comm comm) 
	// throws JpetraException
    {
        this.numGlobalElements = numGlobalElements;
        this.numVnodeElements = numVnodeElements;
        this.elementSize = elementSize;
        this.indexBase = indexBase;
        this.comm = comm;
        minVnodeElementSize = elementSize;
        maxVnodeElementSize = elementSize;
        minElementSize = elementSize;
        maxElementSize = elementSize;
        hasConstantElementSize = true;
        isLinearMap = true;
        isDistributedGlobal = true;
        
        // if (this.numGlobalElements < -1) throw new JpetraException("numGlobalElements < -1");
        // if (this.numVnodeElements < 0) throw new JpetraException("numVnodeElements < 0");
        // if (this.elementSize <= 0) throw new JpetraException("elementSize <= 0");
	if(this.numGlobalElements < -1) {
	    System.out.println("BlockMap constructor: numGlobaElemetns < -1");
	    System.exit(1);
	}
	if(this.numVnodeElements < 0) {
	    System.out.println("BlockMap constructor: numVnodeElements < 0");
	    System.exit(1);
	}
	if(this.elementSize <= 0) {
	    System.out.println("BlockMap constructor: elementSize <= 0");
	    System.exit(1);
	}
        
        int numProc = comm.getNumVnodes();
        int myPID= comm.getVnodeID();
        if (numProc == 1 || numGlobalElements == numVnodeElements) {
            this.numGlobalElements = this.numVnodeElements;
            // Check to see if user's value for numGlobalElements is either -1 
            // (in which case we use our computed value) or matches ours.
            // if (numGlobalElements!=-1 && numGlobalElements!=this.numGlobalElements)
            //    throw new JpetraException("numGlobalElements = " + numGlobalElements
	    //    + ". Should = " + this.numGlobalElements + " or should be set to -1 for me to compute it.");
	    if(numGlobalElements!=-1 && numGlobalElements!=this.numGlobalElements) {
		System.out.println("numGlobalElemetns="+numGlobalElements+". Should="+this.numGlobalElements+" or should be -1 for me to compute it.");
		System.exit(1);
	    }

            if (numProc==1) isDistributedGlobal = false;
            numGlobalEquations = this.numGlobalElements * this.elementSize;
            numVnodeEquations = this.numVnodeElements * this.elementSize;
      
            minAllGID = this.indexBase;
            maxAllGID = minAllGID + this.numGlobalElements - 1;
            minVnodeGID = this.indexBase;
            
            maxVnodeGID = minVnodeGID + this.numVnodeElements - 1;
            minLID = 0;
            maxLID = minLID + this.numVnodeElements - 1;
        }
            
        else if (numProc > 1) {
            // Sum up all local element counts to get global count
            int [] arrayArg0 = new int [1];
            arrayArg0[0] = this.numVnodeElements;
            int [] arrayArg1 = new int[1];
            arrayArg1[0] = this.numGlobalElements;
            arrayArg1 = this.comm.sumAll(arrayArg0);
            this.numGlobalElements = arrayArg1[0];
            
            // Check to see if user's value for numGlobalElements is either -1
            // (in which case we use our computed value) or matches ours.
            // if (numGlobalElements != -1 && numGlobalElements != this.numGlobalElements)
            //    throw new JpetraException("numGlobalElements = " + numGlobalElements
            //    + ". Should = " + this.numGlobalElements + ".");
	    if(numGlobalElements != -1 && numGlobalElements != this.numGlobalElements) {
		System.out.println("numGlobalElements="+numGlobalElements+". Should="+this.numGlobalElements);
		System.exit(1);
	    }
            numGlobalEquations = this.numGlobalElements * this.elementSize;
            numVnodeEquations = this.numVnodeElements* this.elementSize;
            
            minAllGID = this.indexBase;
            maxAllGID = minAllGID + this.numGlobalElements - 1;
            minLID = 0;
            maxLID = minLID + this.numVnodeElements - 1;
            
            // Use the scanSum function to compute a prefix sum of the number of equations
            arrayArg0[0] = this.numVnodeElements;
            arrayArg1[0] = maxVnodeGID;
            arrayArg1 = this.comm.scanSums(arrayArg0);
            maxVnodeGID = arrayArg1[0];
            
            int startIndex = maxVnodeGID - this.numVnodeElements;
            minVnodeGID = startIndex + this.indexBase;
            maxVnodeGID = minVnodeGID + this.numVnodeElements - 1;
        }
            
        //else throw new JpetraException("numGlobalElements = " + numGlobalElements
        //    + ". Should = " + numVnodeElements + " and numProc = " + numProc
        //    + ". Should = 1.");
        if(numGlobalElements != -1 && numGlobalElements != this.numGlobalElements) {
            System.out.println("numGlobalElements="+numGlobalElements+". Should="+this.numGlobalElements+
            " and numProc="+numProc+" should=1");
            System.exit(1);
        }
        
        
        // Setup any information for making global index to local index translation fast.
        globalToLocalSetup();
    }
    
    /**
     * Constructor for a user-defined linear distribution of constant block size elements.
     * @param numGlobalElements the number of elements to distribute
     * @param numVnodeElements     the number of elements owned by the calling process
     * @param globalElements  the global index value of each element on this process;
     *                          length of numVnodeElements not required to be contiguous or
     *                          to be within the range of 0 to numGlobalElements;
     * @param elementSize       the number of equations or vector entries per element
     * @param indexBase         the minimum index value used for arrays that use this
     *                          map, typically 0 for C/C++ and 1 for Fortran
     * @param comm              the communicator containing information on the number
     *                          of processes
     */
    public BlockMap(int numGlobalElements, int numVnodeElements, int [] globalElements,
    int elementSize, int indexBase, Comm comm) 
	//throws JpetraException 
    {
        this.numGlobalElements = numGlobalElements;
        this.numVnodeElements = numVnodeElements;
        this.globalElements = globalElements;
        this.elementSize = elementSize;
        this.indexBase = indexBase;
        this.comm = comm;
        directory = null;
        minVnodeElementSize = elementSize;
        maxVnodeElementSize = elementSize;
	minElementSize = elementSize;
        maxElementSize = elementSize;
        hasConstantElementSize = true;
        isLinearMap = false;
        isDistributedGlobal = true;
        
        int i;
        // Each process gets numVnodeElements equations
        
        // if(this.numGlobalElements < -1) throw new JpetraException("numGlobalElements < -1");
        // if(this.numVnodeElements < 0) throw new JpetraException("numVnodeElements < 0");
        // if(this.elementSize <= 0) throw new JpetraException("elementSize <= 0");        
	if(this.numGlobalElements < -1) {
	    System.out.println("BlockMap constructor: numGlobalElements < -1");
	    System.exit(1);
	}
        if(this.numVnodeElements < 0) {
	    System.out.println("BlockMap constructor: numVnodeElements < 0");
	    System.exit(1);
	}
        if(this.elementSize <= 0) {
	    System.out.println("BlockMap constructor: elementSize <= 0");
	    System.exit(1);
	}
	   
        // Allocate strage for global index list information
        if(numVnodeElements > 0) this.globalElements = new int [numVnodeElements];
        
        // Get procssor information
        int numProc = comm.getNumVnodes();
        int myPID = comm.getVnodeID();
        if(numVnodeElements > 0) {
            // Compute min/max GID on this process
            minVnodeGID = globalElements[0];
            maxVnodeGID = globalElements[0];
            System.arraycopy(globalElements, 0, this.globalElements, 0, numVnodeElements);
            for (i=0; i<numVnodeElements; i++) {
                minVnodeGID = Math.min(minVnodeGID, globalElements[i]);
                maxVnodeGID = Math.max(maxVnodeGID, globalElements[i]);
            }
        }
        
        else {
            minVnodeGID = this.indexBase;
            maxVnodeGID = this.indexBase;
        }
        
        // Local Map and uniprocess case:  Each process gets a complete copy of all elements
        if(numGlobalElements == numVnodeElements || numProc == 1) {
            this.numGlobalElements = this.numVnodeElements;
            // Check to see if user's value for numGlobalElements is either -1 
            // (in which case we use our computed value) or matches ours.
            
            //if(numGlobalElements != -1 && numGlobalElements != this.numGlobalElements)
            //   throw new JpetraException("numGlobalElements = " + numGlobalElements
	    //        + ". Should = " + this.numGlobalElements + " or should be set to -1 for me to compute it.");
            if(numGlobalElements != -1 && numGlobalElements != this.numGlobalElements) {
                System.out.println("numGlobalElements="+numGlobalElements+". Should="+this.numGlobalElements+
                " or should be -1 for me to compute it.");
                System.exit(1);
            }
            
            if(numProc == 1) isDistributedGlobal = false;
            numGlobalEquations = this.numGlobalElements * this.elementSize;
            numVnodeEquations = this.numVnodeElements * this.elementSize;
            
            minAllGID = minVnodeGID;
            maxAllGID = maxVnodeGID;
            minLID = 0;
            maxLID = minLID + this.numVnodeElements - 1;
        }
        
        else if(numProc > 1) {
            // Sum up all local element counts to get global count
            int [] arrayArg0 = new int [1];
            arrayArg0[0] = this.numVnodeElements;
            int [] arrayArg1 = new int [1];
            arrayArg1[0] = this.numGlobalElements;
            arrayArg1 = this.comm.sumAll(arrayArg0);
            this.numGlobalElements = arrayArg1[0];
            // Check to see if user's value for numGlobalElements is either -1
            // (in which case we use our computed value) or matches ours.
            //if(numGlobalElements != -1 && numGlobalElements != this.numGlobalElements)
            //    throw new JpetraException("numGlobalElements = " + numGlobalElements
	    //        + ". Should = " + this.numGlobalElements + " or should be set to -1 for me to compute it.");
            if(numGlobalElements != -1 && numGlobalElements != this.numGlobalElements) {
                System.out.println("numGlobalElements="+numGlobalElements+". Should="+this.numGlobalElements+
                " or should be -1 for me to compute it.");
                System.exit(1);
            }
            
            numGlobalEquations = this.numGlobalElements * this.elementSize;
            numVnodeEquations = this.numVnodeElements * this.elementSize;
            
            minLID = 0;
            maxLID = Math.max(minLID + this.numVnodeElements - 1, minLID);
            
            // Use the Allreduce function to find min/max GID 
            int [] tmpSend = new int [2];
            int [] tmpRecv = new int [2];
            tmpSend[0] = -minVnodeGID; // Negative sign lets us do one reduction
            tmpSend[1] =  maxVnodeGID;
            tmpRecv = this.comm.maxAll(tmpSend);
            minAllGID = -tmpRecv[0];
            maxAllGID =  tmpRecv[1];
            
            //if(minAllGID < this.indexBase) throw new JpetraException(
            //    "minimum element index = " + minAllGID + 
            //    " which is less than the index base = " + this.indexBase + ".");
            if(minAllGID < this.indexBase) {
                System.out.println("minimum element index="+minAllGID+
                " which is less than the index base="+this.indexBase);
                System.exit(1);
            }
        }
            
        //else throw new JpetraException("numGlobalElements = " +
        //    this.numGlobalElements + ". Should = " + this.numVnodeElements +
        //    " and numProc = " + numProc + " should = 1.");
        if(numGlobalElements != -1 && numGlobalElements != this.numGlobalElements) {
            System.out.println("numGlobalElements="+numGlobalElements+". Should="+this.numGlobalElements+
            " and numProc="+numProc+" should=1");
            System.exit(1);
        }

        globalToLocalSetup(); // Setup any information for making global index to local index translation fast.
    }
    
    /**
     * Constructor for a user-defined linear distribution of constant block size elements.
     * @param numGlobalElements the number of elements to distribute
     * @param numVnodeElements     the number of elements owned by the calling process
     * @param globalElements  the global index value of each element on this process;
     *                          length of numVnodeElements not required to be contiguous or
     *                          to be within the range of 0 to numGlobalElements;
     * @param elementSizeList   the element sizes for elements owned by the calling
     *                          process; the ith entry contains the element size of
     *                          the ith element on this procssor
     * @param indexBase         the minimum index value used for arrays that use this
     *                          map, typically 0 for C/C++ and 1 for Fortran
     * @param comm              the communicator containing information on the number
     *                          of processes
     */
    public BlockMap(int numGlobalElements, int numVnodeElements, int [] globalElements,
    int [] elementSizeList, int indexBase, Comm comm) 
	// throws JpetraException 
    {
        
        /* numGlobalElements, numVnodeElements, globalElements, elementSizeList,
         *  indexBase, comm*/
        
        this.numGlobalElements = numGlobalElements;
        this.numVnodeElements = numVnodeElements;
        this.indexBase = indexBase;
        this.comm = comm;
        hasConstantElementSize = false;
        isLinearMap = false;
        isDistributedGlobal = false;
        int i;
        
        // Each process gets numVnodeElements equations
        // if(this.numGlobalElements < -1) throw new JpetraException("numGlobalElements < -1");
        // if(this.numVnodeElements < 0) throw new JpetraException("numVnodeElements < 0");
        // if(this.elementSize <= 0) throw new JpetraException("elementSize <= 0");
        if(this.numGlobalElements < -1) {
	    System.out.println("BlockMap constructor: numGlobalElements < -1");
	    System.exit(1);
	}
        if(this.numVnodeElements < 0) {
	    System.out.println("BlockMap constructor: numVnodeElements < 0");
	    System.exit(1);
	}
	for(i=0; i<this.numVnodeElements; i++) {
	    if(elementSizeList[i] <= 0) {
		System.out.println("BlockMap constructor: elementSize <= 0");
		System.exit(1);
	    }
	}

        for(i=0; i<this.numVnodeElements; i++)
            // if(elementSizeList[i] <= 0) throw new JpetraException("elementSizeList[" + i + "] <= 0.");
	    if(elementSizeList[i] <= 0) {
		System.out.println("BlockMap constructor: elementSizeList["+i+"] <= 0");
		System.exit(1);
	    }
        // Allocate storage for global index list and element size information

        if(numVnodeElements>0) {
            this.globalElements = new int[numVnodeElements];
            this.elementSizeList = new int[numVnodeElements];
        }
        // Get process information

        int numProc = comm.getNumVnodes();
        int myPID = comm.getVnodeID();

        if(numVnodeElements>0) {
            // Compute min/max GID and element size, number of equations on this process
            minVnodeGID = globalElements[0];
            maxVnodeGID = globalElements[0];
            minVnodeElementSize = elementSizeList[0];
            maxVnodeElementSize = elementSizeList[0];
            numVnodeEquations = 0;
            System.arraycopy(globalElements, 0, this.globalElements, 0, numVnodeElements);
            System.arraycopy(elementSizeList, 0, this.elementSizeList, 0, numVnodeElements);
            for(i = 0; i < numVnodeElements; i++) {
                minVnodeGID = Math.min(minVnodeGID,globalElements[i]);
                maxVnodeGID = Math.max(maxVnodeGID,globalElements[i]);
                minVnodeElementSize = Math.min(minVnodeElementSize,elementSizeList[i]);
                maxVnodeElementSize = Math.max(maxVnodeElementSize,elementSizeList[i]);
                numVnodeEquations += elementSizeList[i];
            }
        }
        else {
            minVnodeGID = this.indexBase;
            maxVnodeGID = this.indexBase;
            minVnodeElementSize = 1;
            maxVnodeElementSize = 1;
            this.numVnodeEquations = 0;
        }

        // Local Map and uniprocess case:  Each process gets a complete copy of all elements
        if (numGlobalElements==numVnodeElements || numProc==1) {
            this.numGlobalElements = this.numVnodeElements;
            // Check to see if user's value for numGlobalElements is either -1 
            // (in which case we use our computed value) or matches ours.
            
	    // if (numGlobalElements!=-1 && numGlobalElements!=this.numGlobalElements)
            //    throw new JpetraException("numGlobalElements = " + numGlobalElements
            //    + ". Should = " + this.numGlobalElements + " or should be set to -1 for me to compute it.");
            
            if (numGlobalElements!=-1 && numGlobalElements!=this.numGlobalElements) {
		System.out.println("BlockMap constructor: numGlobalElements="+numGlobalElements+". Should="+this.numGlobalElements+" or be -1 for me to compute it.");
		System.exit(1);
	    }

            if (numProc==1) isDistributedGlobal = false;
            numGlobalEquations = this.numVnodeEquations;

            minAllGID = minVnodeGID;
            maxAllGID = maxVnodeGID;
            minLID = 0;
            maxLID = minLID + this.numVnodeElements - 1;
            minElementSize = minVnodeElementSize;
            maxElementSize = maxVnodeElementSize;
        }
        else if (numProc > 1) {
            // Sum up all local element and equation counts to get global counts
            int [] tmpSend = new int [4];
            int [] tmpRecv = new int [4];
            tmpSend[0] = this.numVnodeElements;
            tmpSend[1] = numVnodeEquations;
            tmpRecv = this.comm.sumAll(tmpSend);
            this.numGlobalElements =  tmpRecv[0];
            numGlobalEquations = tmpRecv[1];

            // Check to see if user's value for numGlobalElements is either -1 
            // (in which case we use our computed value) or matches ours.
            // if (numGlobalElements!=-1 && numGlobalElements!=this.numGlobalElements)
            //    throw new JpetraException("numGlobalElements = " + numGlobalElements
            //    + ". Should = " + this.numGlobalElements + " or should be set to -1 for me to compute it.");
	    if (numGlobalElements!=-1 && numGlobalElements!=this.numGlobalElements) {
		System.out.println("BlockMap constructor: numGlobalElements="+numGlobalElements+". Should="+this.numGlobalElements+" or be -1 for me to compute it.");
		System.exit(1);
	    }

            minLID = 0;
            maxLID = Math.max(minLID + this.numVnodeElements - 1, minLID);
            
            tmpSend[0] = - minVnodeGID; // Negative signs lets us do one reduction
            tmpSend[1] =   maxVnodeGID;
            tmpSend[2] = - minVnodeElementSize;
            tmpSend[3] =   maxVnodeElementSize;

            tmpRecv = this.comm.maxAll(tmpSend);

            minAllGID =      - tmpRecv[0];
            maxAllGID =        tmpRecv[1];
            minElementSize = - tmpRecv[2];
            maxElementSize =   tmpRecv[3];

            // if (minAllGID < this.indexBase) throw new JpetraException(
            //    "minimum element index = " + minAllGID + 
            //    " which is less than the index base = " + this.indexBase + ".");
	    if(minAllGID < this.indexBase) {
		System.out.println("BlockMap constructor: minAllGID is "+minAllGID+" which is less than the index base="+this.indexBase);
		System.exit(1);
	    }
        }
        // else throw new JpetraException("numGlobalElements = " + numGlobalElements
        //    + ". Should = " + numVnodeElements + " and numProc = " + numProc
        //    + ". Should = 1.");
	else {
	    System.out.println("numGlobalElements="+numGlobalElements+". Should="+numVnodeElements+" and numProc="+numProc+". Should=1.");
	    System.exit(1);
	}
      
        globalToLocalSetup(); // Setup any information for making global index to local index translation fast.
    }
    
    public BlockMap(BlockMap map) 
	// throws JpetraException
	{
        this.numGlobalElements = map.numGlobalElements;
        this.numVnodeElements = map.numVnodeElements;
        this.elementSize = map.elementSize;
        this.indexBase = map.indexBase;
        this.comm = map.comm;
        this.numGlobalEquations = map.numGlobalEquations;
        this.numVnodeEquations = map.numVnodeEquations;
        this.minAllGID = map.minAllGID;
        this.maxAllGID = map.maxAllGID;
        this.minVnodeGID = map.minVnodeGID;
        this.maxVnodeGID = map.maxVnodeGID;
        this.minLID = map.minLID;
        this.maxLID = map.maxLID;
        this.minVnodeElementSize = map.minVnodeElementSize;
        this.maxVnodeElementSize = map.maxVnodeElementSize;
        this.minElementSize = map.minElementSize;
        this.maxElementSize = map.maxElementSize;
        this.hasConstantElementSize = map.hasConstantElementSize;
        this.isLinearMap = map.isLinearMap;
        this.isDistributedGlobal = map.isDistributedGlobal;

        int i;
        if (map.globalElements != null) {
            this.globalElements = new int[this.numVnodeElements];
      
            for(i=0; i<this.numVnodeElements; i++)
	        this.globalElements[i] = map.globalElements[i];
        }
        if (map.firstElementEntryList != null) {
            this.firstElementEntryList = new int[this.numVnodeElements+1];
      
            for(i=0; i<numVnodeElements+1; i++)
	        this.firstElementEntryList[i] = map.firstElementEntryList[i];
        }
        if (map.elementSizeList != null) {
            this.elementSizeList = new int[numVnodeElements];
            System.arraycopy(map.elementSizeList, 0, this.elementSizeList, 0, numVnodeElements);
        }
        globalToLocalSetup(); // Setup any information for making global index to local index translation fast.
    }
    
    /**
     * Gets process IDs and corresponding local index value for a given list of
     * global indices.
     *
     * @param numIDs    in on entry; length of GIDList
     * @param GIDList   in on entry; the list of global element numbers
     * @param PIDList   out on exit; the process IDs for GIDList
     * @param LIDList   out on exit; the local index of the GIDList on that process
     */
    public int getRemoteIDList(int numIDs, int [] GIDList, int [] PIDList, 
    int [] LIDList) {
        return getRemoteIDList(numIDs, GIDList, PIDList, LIDList, null);
    }
    
    /**
     * Gets process IDs and corresponding local index value for a given list of
     * global indices.
     *
     * @param NumIDs    in on entry; length of GIDList
     * @param GIDList   in on entry; the list of global element numbers
     * @param PIDList   out on exit; the process IDs for GIDList
     * @param LIDList   out on exit; the local index of the GIDList on that process
     * @param SizeList  out on exit; the element sizes
     */
    public int getRemoteIDList(int numIDs, int [] GIDList, int [] PIDList, 
    int [] LIDList, int [] sizeList) {
        return directory.getDirectoryEntries(numIDs, GIDList, PIDList, LIDList, sizeList);
    }
    
    /**
     * Gets local ID of the passed global ID.
     *
     * @param GID   the global ID
     *
     * @return      the local ID of the passed global ID; -1 if not found on this process
     */
    public int getLID(int GID) {
        if(GID<minVnodeGID || GID > maxVnodeGID) return -1; // Out of range
        else if(!isDistributedGlobal) return (GID-indexBase); // I own all indces
        else if(isLinearMap) return (GID-minVnodeGID); // Can compute with an offset
        else return processLID[GID-minVnodeGID]; // Find it in LID array
    }
    
    /**
     * Gets global ID of passed local ID.
     *
     * @param LID   the local ID
     *
     * @return      the global ID of the passed local ID; indexBase-1 if not found on
     *              this process
     */
    public int getGID(int LID) {
        if(LID < minLID || LID > maxLID) return (indexBase-1); // Out of range
        else if(!isDistributedGlobal) return (LID+indexBase); // I own all indices
        else if(isLinearMap) return (LID+minVnodeGID); // Can compute with an offset
        else return globalElements[LID]; // Find it in globalElements array
    }
    
    /**
     * Gets the LID of the block that contains the given EquationID and the Offset of
     * the equation in that block. ***FIGURE THIS OUT***
     *
     * @param EquationID    the equation of which to get the LID
     * @param BlockID       
     */
    public int findLocalBlockID(int equationID, int [] blockID, int [] blockOffset) {
        int ierr = 0;
        
        if(equationID >= numVnodeEquations) return -1; // Equation is out of range
        
        if(hasConstantElementSize) {
            blockID[0] = equationID / maxElementSize;
            blockOffset[0] = equationID % maxElementSize;
            return 0;
        }
        else {
            if(equationToBlockList == null)
                ierr = getEquationToBlockList(equationToBlockList);
            if(ierr != 0) return ierr;
            if(firstElementEntryList == null)
                ierr = getFirstElementEntryList(firstElementEntryList);
            if(ierr != 0) return ierr;
            
            blockID[0] = equationToBlockList[equationID];
            blockOffset[0] = equationID - firstElementEntryList[blockID[0]];
            return 0;
        }
    }
    
    /**
     * Checks if the passed GID belongs to the calling process in this map.
     *
     * @param GID   the GID to check
     *
     * @return      <code>true</code> if GID belongs to the calling process;
     *              otherwise <code>false</code>
     */
    public boolean isMyGID(int GID) {
        return (getLID(GID) != -1);
    }
    
    /**
     * Checks if the passed LID belongs to the calling process in this map.
     *
     * @param LID   the LID to check
     *
     * @return      <code>true</code> if LID belongs to the calling process;
     *              otherwise <code>false</code>
     */
    public boolean isMyLID(int LID) {
        return (getGID(LID) != indexBase-1);
    }
    
    /**
     * Accessor for the number of elements in the calling process.
     *
     * @return the number of elements in the calling process
     */
    public int getNumVnodeElements() {
        return numVnodeElements;
    }
    
    /**
     * Accessor for the number of global elements.
     *
     * @return the number of global elements
     */
    public int getNumGlobalElements() {
        return numGlobalElements;
    }
    
    /**
     * Accessor for myGlobalElements
     */
    public int getGlobalElements(int[] elements) {
        // If the global element list is not created, then do so.  This can only happen when
        // a linear distribution has been specified.  Thus we can easily construct the update
        // list in this case.
        
        int i;
        if(globalElements == null) {
            for(i=0; i<numVnodeElements; i++)
                elements[i] = minVnodeGID + i;
        }
        else {
            System.arraycopy(globalElements, 0, elements, 0, numVnodeElements);
	}
        
        return 0;
    }
    
    public int[] getGlobalElements() {
        // If elementSizeList not built, do so
        if(globalElements == null && numVnodeElements > 0) {
            int [] tmp = new int [numVnodeElements];
            getGlobalElements(tmp);
            globalElements = tmp;
        }
        return globalElements;
    }

    /**
     * Accessor for the size of each element.
     */
    public int getElementSize() {
        return elementSize;
    }
    
    /**
     * Gets size of element for specified LID
     */
    public int getElementSize(int LID) {
        if(hasConstantElementSize) return elementSize;
        else return(elementSizeList[LID]);
    }
    
    /**
     * Accessor for elementSizeList
     */
    public int getElementSizeList(int [] elementSizeList) {
        // If the element size list is not create, then do so.  This can only happen when
        // a constant element size has been specified.  Thus we can easily construct the block size
        // list in this case.
        
        int i;
        if(this.elementSizeList == null) {
            for(i=0; i<numVnodeElements; i++)
                elementSizeList[i] = elementSize;
	}
        else {
            System.arraycopy(this.elementSizeList, 0, elementSizeList, 0, numVnodeElements);
	}
        return 0;
    }
    
    public int [] getElementSizeList() {
        // If elementSizeList not built, do so
        if(elementSizeList == null && numVnodeElements > 0) {
            int []tmp = new int [numVnodeElements];
            getElementSizeList(tmp);
            elementSizeList = tmp;
        }
        return elementSizeList;
    }
    
    public int getEquationToBlockList(int [] equationToBlockList) {
        // Build an array such that the local block ID is stored for each equation
        
        int i, count = 0;
        if(this.equationToBlockList == null)
            for(i=0; i<numVnodeElements; i++) {
                int size = getElementSize(i);
                for(int j=0; j<size; j++) equationToBlockList[count++] = i;
            }
        else
            System.arraycopy(this.equationToBlockList, 0, equationToBlockList, 0, numVnodeEquations);
        return 0;
    }
    
    public int [] getEquationToBlockList() {
        // If equationToBlockList not built, do so
        if(equationToBlockList == null && numVnodeEquations > 0) {
            int [] tmp = new int [numVnodeEquations];
            getEquationToBlockList(tmp);
            equationToBlockList = tmp;
        }
        return equationToBlockList;
    }
    
    /**
     * Accessor for firstElementEntryList
     */
    public int getFirstElementEntryList(int [] firstElementEntryList) {
        // If the first element entry list is not create, then do so.
        // Note: This array is of length NumMyElement+1
        
        int i;
        
        if(this.firstElementEntryList == null) {
            firstElementEntryList[0] = 0; // First element of first entry is always zero
            
            if(hasConstantElementSize)
                for(i=0; i<numVnodeElements; i++)
                    firstElementEntryList[i+1] = firstElementEntryList[i] + elementSize;
            else
                for(i=0; i<numVnodeElements; i++)
                    firstElementEntryList[i+1] = firstElementEntryList[i] + elementSizeList[i];
        }
        else
            System.arraycopy(this.firstElementEntryList, 0, firstElementEntryList, 0, numVnodeElements);
            
        return 0;
    }
    
    public int [] getFirstElementEntryList() {
        // If elementSizeList not built, do so
        if(firstElementEntryList == null && numVnodeElements > 0) {
            int [] tmp = new int [numVnodeElements+1];
            getFirstElementEntryList(tmp);
            firstElementEntryList = tmp;
        }
        return firstElementEntryList;
    }
    
    /**
     * Accessor for indexBase
     */
    public int getIndexBase() {
        return indexBase;
    }
    
    /**
     * Accessor for numGlobalEquations
     */
    public int getNumGlobalEquations() {
        return numGlobalEquations;
    }
    
    /**
     * Accessor for numMyEquations
     */
    public int getNumVnodeEquations() {
        return numVnodeEquations;
    }
    
    /**
     * Accessor for minAllGID
     */
    public int getMinAllGID() {
        return minAllGID;
    }
    
    /**
     * Accessor for maxAllGID
     */
    public int getMaxAllGID() {
        return maxAllGID;
    }
    
    /**
     * Accessor for minMyGID
     */
    public int getMinVnodeGID() {
        return minVnodeGID;
    }
    
    /**
     * Accessor for maxMyGID
     */
    public int getMaxVnodeGID() {
        return maxVnodeGID;
    }
    
    /**
     * Accessor for minLID
     */
    public int getMinLID() {
        return minLID;
    }
    
    /**
     * Accessor for maxLID
     */
    public int getMaxLID() {
        return maxLID;
    }
    
    /**
     * Accessor for minMyElementSize
     */
    public int getMinVnodeElementSize() {
        return minVnodeElementSize;
    }
    
    /**
     * Accessor for maxMyElementSize
     */
    public int getMaxVnodeElementSize() {
        return maxVnodeElementSize;
    }
    
    /**
     * Accessor for minElementSize
     */
    public int getMinElementSize() {
        return minElementSize;
    }
    
    /**
     * Accessor for maxElementSize
     */
    public int getMaxElementSize() {
        return maxElementSize;
    }
    
    /**
     * Accessor for constantElementSize
     */
    public boolean hasConstantElementSize() {
        return hasConstantElementSize;
    }
    
    /**
     * Compares calling map to passed map
     *
     * @return  <code>true</code> if this and Map are identical maps;
     *          otherwise returns <code>false</code>
     */
    public boolean sameAs(BlockMap Map)
	// throws JpetraException
	{
        if(this == Map) return(true);
    
        if(minAllGID != Map.minAllGID ||
            maxAllGID != Map.maxAllGID ||
            numGlobalElements != Map.numGlobalElements ||
            indexBase != Map.indexBase ) return(false);
  
        if (hasConstantElementSize) {
            if (elementSize != Map.elementSize) return(false);
            else return(true);
        }
        
        else {

            // If we get this far, we need to check local properties and then check across
            // all processes to see if local properties are all true

            int mySameMap = 1; // Assume not needed
            if (numVnodeElements != Map.numVnodeElements) mySameMap = 0;
    
            if (mySameMap==1) 
                for (int i=0; i<numVnodeElements; i++) 
	            if (getGID(i) != Map.getGID(i)) mySameMap = 0;

            // Now get min of MySameMap across all processes

            int globalSameMap = 0;
            int [] arg0 = new int [1];
            arg0[0] = mySameMap;
            int [] arg1 = new int [1];
            arg1[0] = globalSameMap;
            // if(comm.minAll(1, arg0, arg1) != 0) throw new JpetraException
            //    ("Error in function: sameAs");
	    
	    arg1 = comm.minAll(arg0);
	    //if(comm.minAll(1, arg0, arg1) != 0) {
		//System.out.println("BlockMap sameAs: error in call to comm.minAll");
		//System.exit(1);
	    //}
    
            return(globalSameMap == 1);
        }
    }
    
    /**
     * Accessor for linearMap
     */
    public boolean isLinearMap() {
        return isLinearMap;
    }
    
    /**
     * Accessor for distributedGlobal
     */
    public boolean isDistributedGlobal() {
        return isDistributedGlobal;
    }

    /**
     * Accessor for the Comm implementation used by the BlockMap.
     */
    public Comm getComm() {
        return comm;
    }
    
    private void globalToLocalSetup() 
	// throws JpetraException
    {
        int i;
        
        if (numGlobalElements == 0) return; // Nothing to do

        else if (isLinearMap || isDistributedGlobal || numVnodeElements == 0) {
            if (directory == null) directory = new Directory(this); // Make directory
            return; // Nothing else to do
        }
        else {
            // Build processLID vector to make look up of local index values fast
    
            int spanGID = maxVnodeGID - minVnodeGID + 1;
            processLID = new int[spanGID];
    
            for (i=0; i<spanGID; i++) processLID[i] = -1; // Fill all locations with -1
    
            for (i=0; i<numVnodeElements; i++) {
                int tmp = globalElements[i]-minVnodeGID;
                // if (tmp >= 0 || tmp < spanGID) throw new JpetraException("Error in function: globalToLocalSetup");
		if(tmp < 0 || tmp >= spanGID) {
		    System.out.println("Error in function globalToLocalSetup");
		    System.exit(1);
		}

                // assert(tmp>=0); assert(tmp <SpanGID);
                processLID[globalElements[i]-minVnodeGID] = i; // Spread local indices
            }
    
            if (directory == null) directory = new Directory(this); // Make directory
        }
    }
    
    public Object clone() throws CloneNotSupportedException {
        return super.clone();
    } 

    public int [] GETGE() {
	return globalElements;
    }
    public int [] GETESL() {
	return elementSizeList;
    }
}
