// @HEADER
// ***********************************************************************
// 
//               Java Implementation of the Petra Library
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

/*
 * MultiVector.java
 *
 * Created on June 14, 2001, 11:35 AM
 */

package Jpetra;

import org.netlib.util.*;
import org.netlib.blas.*;

/**
 *
 * @author  Michael William Boldt
 * @version 
 */
public class MultiVector extends JpetraObject {
    
    protected int indexBase = 0;
    protected double [][] values = null;
    private BlockMap map;
    protected double flops = 0.0;
    private int vnodeLength = 0;
    
    // User-queriable attributes
    private int globalLength = 0;
    private int numVectors = 0;
    private boolean userAllocated = false;
    private boolean distributedGlobal = false;
    private boolean allocated = false;
    private double seed = 1.0;
    private double [] imports = null;
    private double [] exports = null;
    private int lenImporters = 0;
    private int lenExporters = 0;
    
    // Global reduction variables
    private Comm comm;
    private double [] doubleTemp = null;
    private int [] intTemp = null;

    /**
     * Creates a <code>MultiVector</code> object and fills it with zero values.
     *
     * @param map           A <code>LocalMap</code>, <code>Map</code>, or 
                            <code>BlockMap</code>
     * @param numVectors    Number of vectors in the multi-vector
     */
    public MultiVector(BlockMap map, int numVectors) /*throws JpetraException*/ {
        
        initialize(map, numVectors);
        allocate("Copy");
        putScalar(0.0); // Fill all vectors with zero
    }
    
    /**
     * Creates a copy of the passed <code>MultiVector</code>.
     *
     * @param source        The multi-vector of which to make a copy
     */
    public MultiVector(MultiVector source) /*throws JpetraException*/ {
        
        initialize(source.map, source.numVectors);
        allocate("Copy");
        doCopy(source.values);
    }
    
    /**
     * Creates a <code>MultiVector</code> and sets the values from the
     * gives it the values from <code>values</code>.
     *
     * @param copyView      "Copy" makes a deep copy of the values;
     *                      "View" makes a shallow copy of the values;
     *                      Any other string thows an exception
     * @param map           A <code>LocalMap</code>, <code>Map</code>, or 
     *                      <code>BlockMap</code>
     * @param values        The values to be copied
     * @param numVectors    Number of vectors in multi-vector
     */
    public MultiVector(String copyView, BlockMap map, double [][] values,
		       int numVectors) /*throws JpetraException*/ {
        
        initialize(map, numVectors);
        allocate(copyView);
        if(copyView == "Copy") doCopy(values);
        else doView(values);
    }
    
    public MultiVector(String copyView, MultiVector source, 
		       int [] indices, int numVectors) /*throws JpetraException*/ {
        
        initialize(source, numVectors);
        allocate(copyView);
        
        double [][] sourceValues = new double[numVectors][];
        for(int i=0; i<numVectors; i++) sourceValues[i] = source.values[indices[i]];
        
        if(copyView == "Copy") doCopy(sourceValues);
        else doView(sourceValues);
    }
    
    /**
     * Set multi-vector values from range of vectors in an existing <code>
     * MultiVector</code>.
     *
     * @param copyView      "Copy" makes a deep copy of the values;
     *                      "View" makes a shallow copy of the values;
     *                      Any other string throws an exception
     * @param source        An existing fully constructed <code>MultiVector</code>
     * @param startIndex    First of the vectors to copy
     * @param numVectors    Number of vectors in new multi-vector
     */
    public MultiVector(String copyView, MultiVector source,
		       int startIndex, int numVectors) /*throws JpetraException*/ {
        
        initialize(source, numVectors);
        allocate(copyView);
        
        double[][] sourceValues = new double [numVectors] [vnodeLength];
        System.arraycopy(source, startIndex, sourceValues, 0, numVectors);
        if(copyView == "Copy") doCopy(sourceValues);
        else doView(sourceValues);
    }
    
    /**
     * Initializes many member variables; Called from constructor
     */
    private void initialize(BlockMap map, int numVectors) {
        indexBase = map.getIndexBase();
        this.map = map;
        vnodeLength = map.getNumMyEquations();
        globalLength = map.getNumGlobalEquations();
        this.numVectors = numVectors;
        distributedGlobal = map.isDistributedGlobal();
        comm = map.getComm();
    }
    
    /**
     * Initializes many member variables; Called from constructor
     */
    private void initialize(MultiVector source, int numVectors) {
        indexBase = source.indexBase;
        map = source.map;
        vnodeLength = source.vnodeLength;
        globalLength = source.globalLength;
        this.numVectors = numVectors;
        distributedGlobal = source.distributedGlobal;
        comm = source.comm;
    }
    
    /**
     * Allocates proper space for <code>values</code>, <code>doubleTemp</code>,
     * and <code>intTemp</code>.
     *
     * @param copyView      "Copy" allocates space for a deep copy;
     *                      "View" allocates space for a shallow copy;
     *                      any other string thows an exception
     */
    private int allocate(String copyView) /*throws JpetraException*/ {
        if(allocated) return 0;
        
	/*
        if(numVectors <= 0) throw new JpetraException
            ("number of vectors must be greater than zero, is " + numVectors);
	*/
	if(numVectors <= 0) {
	    System.out.println("Number of vectors must be greater than 0, is "+numVectors);
	    System.exit(1);
	}
        if(copyView == "Copy") values = new double [numVectors][vnodeLength];
        else if(copyView == "View") values = new double [numVectors][];
	/*
        else throw new JpetraException
            ("argument must be either \"Copy\" or \"View\"");
	*/
	else {
	    System.out.println("Argument must be either \"Copy\" or \"View\"");
	    System.exit(1);
	}
        
        doubleTemp = new double [numVectors];
        intTemp = new int [numVectors];
        if(distributedGlobal) seed = 2 * comm.getVnodeID() + 1;
        
        allocated = true;
        userAllocated = true;
        
        return 0;
    }
    
    /**
     * Puts a deep copy of <code>source</code> in <code>values</code>.
     *
     * @param source        The values to copy
     */
    private void doCopy(double [][] source) {
        for(int i=0; i<numVectors; i++)
            System.arraycopy(source[i], 0, values[i], 0, vnodeLength);
    }
    
    /**
     * Puts a shallow copy of <code>source</code> in <code>values</code>.
     *
     * @param source        The values to copy
     */
    private void doView(double [][] source) {
	System.arraycopy(source, 0, values, 0, numVectors);
    }
    
    /**
     * Generates random numbers with an approximately uniform distribution
     * in the interval (-1,1).
     */
    public int random() {
        // For MPI change so not all processes get same randoms
        int i, j;
        for(i=0; i<numVectors; i++)
            for(j=0; j<vnodeLength; j++)
                values[i][j] = 2 * Math.random() - 1;
        return 0;
    }
    
    /**
     * Extracts a copy of a <code>MultiVector</code> into a user's
     * two-dimensional array.
     */
    public int extractCopy(double [][] vectorArray) {
        for(int i=0; i<numVectors; i++)
            System.arraycopy(values[i], 0, vectorArray[i], 0, vnodeLength);
        return 0;
    }

    public int extractView(double [][] vectorArray) {
	vectorArray = values;
	return 0;
    }

    public double [][] extractView() {
	return values;
    }
    
    /**
     * Fills multi-vector with the value <code>scalar</code>
     */
    public int putScalar(double scalar) {
        for(int i=0; i<numVectors; i++)
            for(int j=0; j<vnodeLength; j++) values[i][j] = scalar;
        return 0;
    }
    
    /**
     * Importers a set of values from the input multi-vector.
     *
     * @param a             Multi-vector from which values are imported to calling multi-vector
     *
     * @param importer      Specifies the communication required
     *
     * @param combineMode   Specivies how results should be combined on the receiving
     *                      process.
     */
    public int doImporter(MultiVector a, Importer importer, String combineMode)
    throws JpetraException {
        if(numVectors != a.numVectors) return -1;
        if(!map.sameAs(importer.getTargetMap())) return -2;
        if(values != a.values && !a.map.sameAs(importer.getSourceMap())) return -3;
        
        int numSameIDs = importer.getNumSameIDs();
        int numPermuteIDs = importer.getNumPermuteIDs();
        int numRemoteIDs = importer.getNumRemoteIDs();
        int numExporterIDs = importer.getNumExporterIDs();
        int [] exportLIDs = importer.getExporterLIDs();
        int [] remoteLIDs = importer.getRemoteLIDs();
        int [] permuteToLIDs = importer.getPermuteToLIDs();
        int [] permuteFromLIDs = importer.getPermuteFromLIDs();
        int nSend = numVectors * importer.getNumSend();
        int nRecv = numVectors * importer.getNumRecv();
        
        return doTransfer(a, combineMode, numSameIDs, numPermuteIDs, numRemoteIDs,
        numExporterIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
        nSend, nRecv, lenExporters, exports, lenImporters, imports, false);
    }
    
    /**
     * Importers a set of values from the input multi-vector using the <code>Importer</code>
     * object.
     *
     * @param a             Multi-vector from which values are imported to calling multi-vector
     *
     * @param exporter      Specifies the communication required
     *
     * @param combineMode   Specivies how results should be combined on the receiving
     *                      process.
     */
    public int doImporter(MultiVector a, Exporter exporter, String combineMode)
    throws JpetraException {
        if(numVectors != a.numVectors) return -1;
        if(!map.sameAs(exporter.getTargetMap())) return -2;
        if(values != a.values && !a.map.sameAs(exporter.getSourceMap())) return -3;
        
        int numSameIDs = exporter.getNumSameIDs();
        int numPermuteIDs = exporter.getNumPermuteIDs();
        int numRemoteIDs = exporter.getNumRemoteIDs();
        int numExporterIDs = exporter.getNumExporterIDs();
        int [] exportLIDs = exporter.getExporterLIDs();
        int [] remoteLIDs = exporter.getRemoteLIDs();
        int [] permuteToLIDs = exporter.getPermuteToLIDs();
        int [] permuteFromLIDs = exporter.getPermuteFromLIDs();
        int nSend = numVectors * exporter.getNumSend();
        int nRecv = numVectors * exporter.getNumRecv();
        
        return doTransfer(a, combineMode, numSameIDs, numPermuteIDs, numRemoteIDs,
        numExporterIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
        nSend, nRecv, lenExporters, exports, lenImporters, imports, false);
    }
    
    /**
     * Exporters a set of values from the input multi-vector.
     *
     * @param a             Multi-vector from which values are imported to calling multi-vector
     *
     * @param importer      Specifies the communication required
     *
     * @param combineMode   Specivies how results should be combined on the receiving
     *                      process.
     */
    public int doExporter(MultiVector a, Importer importer, String combineMode)
    throws JpetraException {
        if(numVectors != a.numVectors) return -1;
        if(!map.sameAs(importer.getTargetMap())) return -2;
        if(values != a.values && !a.map.sameAs(importer.getSourceMap())) return -3;
        
        int numSameIDs = importer.getNumSameIDs();
        int numPermuteIDs = importer.getNumPermuteIDs();
        int numRemoteIDs = importer.getNumRemoteIDs();
        int numExporterIDs = importer.getNumExporterIDs();
        int [] exportLIDs = importer.getExporterLIDs();
        int [] remoteLIDs = importer.getRemoteLIDs();
        int [] permuteToLIDs = importer.getPermuteToLIDs();
        int [] permuteFromLIDs = importer.getPermuteFromLIDs();
        int nSend = numVectors * importer.getNumSend();
        int nRecv = numVectors * importer.getNumRecv();
        
        return doTransfer(a, combineMode, numSameIDs, numPermuteIDs, numRemoteIDs,
        numExporterIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
        nSend, nRecv, lenImporters, imports, lenExporters, exports, true);
    }
    
    /**
     * Exporters a set of values from the input multi-vector.
     *
     * @param a             Multi-vector from which values are imported to calling multi-vector
     *
     * @param exporter      Specifies the communication required
     *
     * @param combineMode   Specivies how results should be combined on the receiving
     *                      process.
     */
    public int doExporter(MultiVector a, Exporter exporter, String combineMode)
    throws JpetraException {
        if(numVectors != a.numVectors) return -1;
        if(!map.sameAs(exporter.getTargetMap())) return -2;
        if(values != a.values && !a.map.sameAs(exporter.getSourceMap())) return -3;
        
        int numSameIDs = exporter.getNumSameIDs();
        int numPermuteIDs = exporter.getNumPermuteIDs();
        int numRemoteIDs = exporter.getNumRemoteIDs();
        int numExporterIDs = exporter.getNumExporterIDs();
        int [] exportLIDs = exporter.getExporterLIDs();
        int [] remoteLIDs = exporter.getRemoteLIDs();
        int [] permuteToLIDs = exporter.getPermuteToLIDs();
        int [] permuteFromLIDs = exporter.getPermuteFromLIDs();
        int nSend = numVectors * exporter.getNumSend();
        int nRecv = numVectors * exporter.getNumRecv();
        
        return doTransfer(a, combineMode, numSameIDs, numPermuteIDs, numRemoteIDs,
        numExporterIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
        nSend, nRecv, lenImporters, imports, lenExporters, exports, true);
    }
    
    private int doTransfer(MultiVector a, String combineMode, int numSameIDs,
    int numPermuteIDs, int numRemoteIDs, int numExporterIDs, int [] permuteToLIDs,
    int [] permuteFromLIDs, int [] remoteLIDs, int [] exportLIDs, int nSend,
    int nRecv, int lenExporters, double [] exports, int lenImporters,
    double [] imports, boolean doReverse) {
        int ierr = 0;
        double [][] aValues = a.values;
        
        int [] firstElementEntryList = null;
        int [] aFirstElementEntryList = null;
        int [] elementSizeList = null;
        int [] aElementSizeList = null;
        if(!map.hasConstantElementSize()) {
            firstElementEntryList = map.getFirstElementEntryList();
            elementSizeList = map.getElementSizeList();
        }
        if(!a.map.hasConstantElementSize()) {
            aFirstElementEntryList = a.map.getFirstElementEntryList();
            aElementSizeList = a.map.getElementSizeList();;
        }
        ierr = copyAndPermute(values, aValues, numVectors, a.map.getMaxElementSize(),
        a.map.hasConstantElementSize(), aElementSizeList, aFirstElementEntryList,
        firstElementEntryList, numSameIDs, numPermuteIDs, permuteToLIDs, permuteFromLIDs);
        
        if(ierr != 0) return ierr;
        if(combineMode == "Zero") return 0;
        
        if(nSend > lenExporters) {
            exports = new double [nSend];
            lenExporters = nSend;
        }
        
        ierr = pack(aValues, numVectors, a.map.getMaxElementSize(),
        a.map.hasConstantElementSize(), aElementSizeList, aFirstElementEntryList,
        numExporterIDs, exportLIDs, exports);
        
        if(ierr != 0) return ierr;
        
        if(nRecv > lenExporters) {
            imports = new double [nRecv];
            lenImporters = nRecv;
        }
        
        // #ifdef PETRA_MPI ...
        
        return ierr;
    }
    
    private int copyAndPermute(double [][] to, double [][] from, int numVectors,
    int maxElementSize, boolean constantElementSize, int [] fromElementSizeList,
    int [] fromFirstElementEntryList, int [] toFirstElementEntryList, int numSameIDs,
    int numPermuteIDs, int [] permuteToLIDs, int [] permuteFromLIDs) {
        
        int i, j, jj, jjj, k;
        int numSameEntries;
        boolean case1 = false;
        boolean case2 = false;
        
        if(maxElementSize == 1) {
            case1 = true;
            numSameEntries = numSameIDs;
        }
        else if(constantElementSize) {
            case2 = true;
            numSameEntries = numSameIDs * maxElementSize;
        }
        else numSameEntries = fromFirstElementEntryList[numSameIDs];
        
        // Short circuit for the case where the source and target vector is the same
        if(to == from) numSameEntries = 0;
        
        // Do copy first
        if(numSameIDs > 0)
            if(to != from)
                for(i=0; i<numVectors; i++)
                    for(j=0; j<numSameEntries; j++)
                        to[i][j] = from[i][j];
        
        // Do local permutation next
        if(numPermuteIDs > 0)
            // Point entry case
            if(case1)
                if(numVectors == 1)
                    for(j=0; j<numPermuteIDs; j++)
                        to[0][permuteToLIDs[j]] = from[0][permuteFromLIDs[j]];
                else
                    for(j=0; j<numPermuteIDs; j++) {
                        jj = permuteToLIDs[j];
                        jjj = permuteFromLIDs[j];
                        for(i=0; i<numVectors; i++)
                            to[i][jj] = from[i][jjj];
                    }
            else if(case2)
                for(j=0; j<numPermuteIDs; j++) {
                    jj = maxElementSize * permuteToLIDs[j];
                    jjj = maxElementSize * permuteFromLIDs[j];
                    for(i=0; i<numVectors; i++)
                        for(k=0; k<maxElementSize; k++)
                            to[i][jj+k] = from[i][jjj+k];
                }
            else
                for(j=0; j<numPermuteIDs; j++) {
                    jj = toFirstElementEntryList[permuteToLIDs[j]];
                    jjj = fromFirstElementEntryList[permuteFromLIDs[j]];
                    int elementSize = fromElementSizeList[permuteFromLIDs[j]];
                    for(i=0; i<numVectors; i++)
                        for(k=0; k<elementSize; k++)
                            to[i][jj+k] = from[i][jjj+k];
                }
        return 0;
    }
    
    private int pack(double [][] from, int numVectors, int maxElementSize,
    boolean constantElementSize, int [] fromElementSizeList, int [] fromFirstElementEntryList,
    int numSendIDs, int [] sendLIDs, double [] sends) {
        int i, j, jj, k, count = 0;
        
        double [] tmp;
        
        if(numSendIDs > 0) {
            tmp = sends;
            
            if(maxElementSize == 1) {
                count = 0;
                if(numVectors == 1) for(j=0; j<numSendIDs; j++)
                    tmp[j] = from[0][sendLIDs[j]];
                else
                    for(j=0; j<numSendIDs; j++) {
                        jj = sendLIDs[j];
                        for(i=0; i<numVectors; i++)
                            tmp[count++] = from[i][jj];
                    }
            }
            
            else if(constantElementSize) {
                count = 0;
                for(j=0; j<numSendIDs; j++) {
                    jj = maxElementSize * sendLIDs[j];
                    for(i=0; i<numVectors; i++)
                        for(k=0; k<maxElementSize; k++)
                            tmp[count++] = from[i][jj+k];
                }
            }
            
            else {
                count = 0;
                int sizeOfPacket = numVectors * maxElementSize;
                for(j=0; j<numSendIDs; j++) {
                    jj = fromFirstElementEntryList[sendLIDs[j]];
                    int elementSize = fromElementSizeList[sendLIDs[j]];
                    for(i=0; i<numVectors; i++)
                        for(k=0; k<elementSize; k++)
                            sends[count++] = from[i][jj+k];
                }
            }
        }
        return 0;
    }
    
    private int unpackAndCombine(double [][] to, int numVectors, int maxElementSize,
    boolean constantElementSize, int [] toElementSizeList, int [] toFirstElementEntryList,
    int numRecvIDs, int [] recvLIDs, double [] recvs, String combineMode) throws JpetraException {
        int i, j, jj, k, p = 0;
        double [] tmp;
        
        // Unpack it
        if(numRecvIDs > 0) {

            // Point entry case
            if(maxElementSize == 1) {
                if(numVectors == 1) {
                    if(combineMode == "Add")
                        for(j=0; j<numRecvIDs; j++) to[0][recvLIDs[j]] += recvs[p++];
                    else if(combineMode == "Insert")
                        for(j=0; j<numRecvIDs; j++) to[0][recvLIDs[j]] = recvs[p++];
                    else if(combineMode == "Average")
                        for(j=0; j<numRecvIDs; j++) {
                            to[0][recvLIDs[j]] += recvs[p++];
                            to[0][recvLIDs[j]] *= 0.5;
                        }
                    else throw new JpetraException
                    ("combineMode must be \"Add\", \"Insert\", or \"Average\"");
                }
                else { // numVectors > 1
                    if(combineMode == "Add") {
                        for(j=0; j<numRecvIDs; j++) {
                            jj = recvLIDs[j];
                            for(i=0; i<numVectors; i++)
                                to[i][j] = recvs[p++];
                        }
                    }
                    else if(combineMode == "Insert") {
                        for(j=0; j<numRecvIDs; j++) {
                            jj = recvLIDs[j];
                            for(i=0; i<numVectors; i++)
                                to[i][jj] = recvs[p++];
                        }
                    }
                    else if(combineMode == "Average") {
                        for(j=0; j<numRecvIDs; j++) {
                            jj = recvLIDs[j];
                            for(i=0; i<numVectors; i++) {
                                to[i][jj] += recvs[p++];
                                to[i][jj] *= 0.5;
                            }
                        }
                    }
                    else throw new JpetraException
                    ("combineMode must be \"Add\", \"Insert\", or \"Average\"");
                }
            }
            
            // Constant block size case
            
            else if(constantElementSize) {
                if(combineMode == "Add") {
                    for(j=0; j<numRecvIDs; j++) {
                        jj = maxElementSize*recvLIDs[j];
                        for(i=0; i<numVectors; i++)
                            for(k=0; k<maxElementSize; k++)
                                to[i][jj+k] += recvs[p++];
                    }
                }
                else if(combineMode == "Insert") {
                    for(j=0; j<numRecvIDs; j++) {
                        jj = maxElementSize * recvLIDs[j];
                        for(i=0; i<numVectors; i++)
                            for(k=0; k<maxElementSize; k++)
                                to[i][jj+k] = recvs[p++];
                    }
                }
                else if(combineMode == "Average") {
                    for(j=0; j<numRecvIDs; j++) {
                        jj = maxElementSize * recvLIDs[j];
                        for(i=0; i<numVectors; i++)
                            for(k=0; k<maxElementSize; k++) {
                                to[i][jj+k] += recvs[p++];
                                to[i][jj+k] *= 0.5;
                            }
                    }
                }
                else throw new JpetraException
                ("combineMode must be \"Add\", \"Insert\", or \"Average\"");
            }
            else {
                int sizeOfPacket = numVectors * maxElementSize;
                if(combineMode == "Add") {
                    for(j=0; j<numRecvIDs; j++) {
                        p = j * sizeOfPacket;
                        jj = toFirstElementEntryList[recvLIDs[j]];
                        int elementSize = toElementSizeList[recvLIDs[j]];
                        for(i=0; i<numVectors; i++)
                            for(k=0; k<elementSize; k++)
                                to[i][jj+k] += recvs[p++];
                        }
                }
                else if(combineMode == "Insert") {
                    for(j=0; j<numRecvIDs; j++) {
                        p = j * sizeOfPacket;
                        jj = toFirstElementEntryList[recvLIDs[j]];
                        int elementSize = toElementSizeList[recvLIDs[j]];
                        for(i=0; i<numVectors; i++)
                            for(k=0; k<elementSize; k++)
                                to[i][jj+k] = recvs[p++];
                    }
                }
                else if(combineMode == "Average") {
                    for(j=0; j<numRecvIDs; j++) {
                        p = j * sizeOfPacket;
                        jj = toFirstElementEntryList[recvLIDs[j]];
                        int elementSize = toElementSizeList[recvLIDs[j]];
                        for(i=0; i<numVectors; i++)
                            for(k=0; k<elementSize; k++) {
                                to[i][jj+k] += recvs[p++];
                                to[i][jj+k] *= 0.5;
                            }
                    }
                }
                else throw new JpetraException
                ("combineMode must be \"Add\", \"Insert\", or \"Average\"");
            }
        }
        return 0;
    }
    
    /*
    public int dot(MultiVector a, double [] result) {
        
        int i;
        
        if(numVectors != a.numVectors) return -1;
        if(vnodeLength != a.vnodeLength) return -2;
        double [] localResult = new double [numVectors];
        
        for(i=0; i<numVectors; i++) {
            sum = 0.0;
            for(j=0; j<vnodeLength; j++) sum += values[i][j] * a.values[i][j];
            localResult[i] = sum;
        }
        
        return map.getComm().sumAll(numVectors, localResult, result);
    }
    */
    public int dot(MultiVector a, double [] result) {
        int i;
        
        if(numVectors != a.numVectors) return -1;
        if(vnodeLength != a.vnodeLength) return -2;
        
        double [][] aValues = a.values;
        
        for(i=0; i<numVectors; i++) {
	    doubleTemp[i] = 0.0;
	    for (int j=0; j<vnodeLength; j++) doubleTemp[i] += values[i][j]*aValues[i][j];
	} 
	    // doubleTemp[i] = DDOT.DDOT(vnodeLength, values[i], 0, aValues[i], 0);
        result = comm.sumAll(doubleTemp);
        
        updateFlops(2 * globalLength * numVectors);
        
        return 0;
    }
    
    /**
     * Puts element-wise absolute values of input multi-vector in calling multi-vector.
     */
    public int abs(MultiVector a) {
        
        int i, j;
        
        if(numVectors != a.numVectors) return -1;
        if(vnodeLength != a.vnodeLength) return -2;
        
        double [][] aValues = a.values;
        
        for(i=0; i<numVectors; i++)
            for(j=0; j<vnodeLength; j++)
                values[i][j] = Math.abs(aValues[i][j]);
        return 0;
    }
    
    /**
     * Puts element-wise reciprocal values of input multi-vector in calling multi-vector.
     */
    public int reciprocal(MultiVector a) {
        int i, j, length;
        double [][] aValues = a.values;
        
        for(i=0; i<numVectors; i++) {
            length = values[i].length;
            for(j=0; j<length; j++)
                if(aValues[i][j] <= Double.MIN_VALUE && aValues[i][j] >= -Double.MIN_VALUE)
                    values[i][j] = Double.MAX_VALUE;
                else values[i][j] = 1.0 / aValues[i][j];
        }
        return 0;
    }
    
    /**
     * Scales values of calling multi-vector by <code>scalar</code>.
     */
    public int scale(double scalar) {
        int i;
        
        for(i=0; i<numVectors; i++)
            DSCAL.DSCAL(vnodeLength, scalar, values[i], 1);
        
        updateFlops(globalLength * numVectors);
        
        return 0;
    }
    
    /**
     * Replaces values of calling multi-vector with scaled values of <code>a</code>.
     */
    public int scale(double scalar, MultiVector a) {
        
        int i, j;
        
        if(numVectors != a.numVectors) return -1;
        if(vnodeLength != a.vnodeLength) return -2;
        
        double [][] aValues = a.values;
        
        for(i=0; i<numVectors; i++)
            for(j=0; j<vnodeLength; j++) values[i][j] = scalar * aValues[i][j];
        updateFlops(globalLength * numVectors);
        return 0;
    }
    
    /**
     * this = scalar*this + scalarA*a
     */
    public int update(double scalarA, MultiVector a, double scalar) {
        int i, j;
        
        if(numVectors != a.numVectors) return -1;
        if(vnodeLength != a.vnodeLength) return -2;
        
        double [][] aValues = a.values;
        
        if(scalar == 0.0) {
            for(i=0; i<numVectors; i++)
                for(j=0; j<vnodeLength; j++)
                    values[i][j] = scalarA * aValues[i][j];
            updateFlops(globalLength * numVectors);
        }
        
        else if(scalar == 1.0) {
            for(i=0; i<numVectors; i++)
                for(j=0; j<vnodeLength; j++)
                    values[i][j] = values[i][j] + scalarA * aValues[i][j];
            updateFlops(2 * globalLength * numVectors);
        }
        else if(scalarA == 1.0) {
	    for(i=0; i<numVectors; i++)
		for(j=0; j<vnodeLength; j++)
		    values[i][j] = scalar * values[i][j] + aValues[i][j];
	    updateFlops(2 * globalLength * numVectors);
	}
        else {
            for(i=0; i<numVectors; i++)
                for(j=0; j<vnodeLength; j++)
                    values[i][j] = scalar * values[i][j] + scalarA * aValues[i][j];
            updateFlops(3 * globalLength * numVectors);
        }
        
        return 0;
    }
    
    /**
     * this = scalar*this + scalarA*a + scalarB*b
     */
    public int update(double scalarA, MultiVector a,
    double scalarB, MultiVector b, double scalar) {
        int i, j;
        
        if(scalarA == 0.0) return update(scalarB, b, scalar);
        if(scalarB == 0.0) return update(scalarA, a, scalar);
        
        if(numVectors != a.numVectors || numVectors != b.numVectors) return -1;
        if(vnodeLength != a.vnodeLength || vnodeLength != b.vnodeLength) return -2;
        
        double [][] aValues = a.values;
        double [][] bValues = b.values;
        
        if(scalar == 0.0) {
            if(scalarA == 1.0) {
                for(i=0; i<numVectors; i++)
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = aValues[i][j] + scalarB * bValues[i][j];
                updateFlops(2 * globalLength * numVectors);
            }
            
            else if(scalarB == 1.0) {
                for(i=0; i<numVectors; i++)
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = bValues[i][j] + scalarA * aValues[i][j];
                updateFlops(2 * globalLength * numVectors);
            }
            
            else {
                for(i=0; i<numVectors; i++)
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = scalarA * aValues[i][j] + scalarB * bValues[i][j];
                updateFlops(3 * globalLength * numVectors);
            }
        }
        
        else if(scalar == 1.0) {
            if(scalarA == 1.0) {
                for(i=0; i<numVectors; i++)
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] += aValues[i][j] + scalarB * bValues[i][j];
                updateFlops(3 * globalLength * numVectors);
            }
            
            else if(scalarB == 1.0) {
                for(i=0; i<numVectors; i++)
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] += bValues[i][j] + scalarA * aValues[i][j];
                updateFlops(3 * globalLength * numVectors);
            }
            
            else {
                for(i=0; i<numVectors; i++)
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] += scalarA * aValues[i][j] + scalarB * bValues[i][j];
                updateFlops(4 * globalLength * numVectors);
            }
        }
        
        else {
            if(scalarA == 1.0) {
                for(i=0; i<numVectors; i++)
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = scalar * values[i][j] + aValues[i][j] + scalarB * bValues[i][j];
                updateFlops(4 * globalLength * numVectors);
            }
            
            else if(scalarB == 1.0) {
                for(i=0; i<numVectors; i++)
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = scalar * values[i][j] + scalarA * aValues[i][j] + bValues[i][j];
                updateFlops(4 * globalLength * numVectors);
            }
            
            else {
                for(i=0; i<numVectors; i++)
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = scalar * values[i][j] + scalarA * aValues[i][j] + scalarB * bValues[i][j];
                updateFlops(5 * globalLength * numVectors);
            }
        }
        
        return 0;
    }
    
    /**
     * Computes the 1-norm of ech vector in calling multi-vector.
     *
     * @param result    result[i] will contain the 1-norm of the ith vactor
     */
    public int norm1(double [] result) {
        if(result == null) result = new double [numVectors];
        else if(result.length != numVectors) return -1;
        
        int i, j;
        double sum;
 
        for(i=0; i<numVectors; i++) {
            //doubleTemp[i] = DASUM.DASUM(vnodeLength, values[i], 1);
            sum = 0.0;
            for(j=0; j<vnodeLength; j++) sum += Math.abs(values[i][j]);
            doubleTemp[i] = sum;
        }

        result = comm.sumAll(doubleTemp);
        
        updateFlops(2 * globalLength * numVectors);
        
        return 0;
    }
    
    /**
     * Computes the 2-norm of each vector in calling multi-vector.
     *
     * @param result    result[i] will contain the 1-norm of the ith vector
     */
    public int norm2(double [] result) {
        if(result == null) result = new double [numVectors];
        else if(result.length != numVectors) return -1;
        
        int i, j;
        double sum;
        
        for(i=0; i<numVectors; i++) {
            sum = 0.0;
            for(j=0; j<vnodeLength; j++) sum += values[i][j] * values[i][j];
            doubleTemp[i] = sum;
        }
        result = comm.sumAll(doubleTemp);
        for(i=0; i<numVectors; i++) result[i] = Math.sqrt(result[i]);
        updateFlops(2 * globalLength * numVectors);
        return 0;
    }
    
    public int normInf(double [] result) {
        if(result == null) result = new double [numVectors];
        else if(result.length != numVectors) return -1;
        
        int i, j;
        double maxval;
 
        for(i=0; i<numVectors; i++) {
            maxval = 0.0;
            for(j=0; j<vnodeLength; j++) {
                double temp = Math.abs(values[i][j]);
                if (temp>maxval) maxval = temp;
            }
            doubleTemp[i] = maxval;
            //j = IDAMAX.IDAMAX(vnodeLength, values[i], 1) - 1;
            //doubleTemp[i] = Math.abs(values[i][j]);
        }
 
        result = comm.maxAll(doubleTemp);
       
        return 0;
    }
    
    /**
     * Computes weighted 2-norm (RMS norm) of each vector in calling multi-vector.
     *
     * @param weights   Multi-vector of weights; If <code>weights</code> contains 
     * a single vector, that vector will be used as the weights for all vectors, 
     * otherwise should have the same number of vectors as calling multi-vector
     *
     * @param result    result[i] will contain the weighted 2-norm of the ith vector
     */
    public int normWeighted(MultiVector weights, double [] result) {

        int i, j;
        double sum;
	boolean oneW = false;

	if(weights.numVectors == 1) oneW = true;
	else if(numVectors != weights.numVectors) return -1;
        if(vnodeLength != weights.vnodeLength) return -2;
        
        double [][] wValues = weights.values;
        
        if(oneW) {
            double [] w = wValues[0];
            for(i=0; i<numVectors; i++) {
                sum = 0.0;
                for(j=0; j<vnodeLength; j++) {
                    double tmp = values[i][j] / w[j];
                    sum += tmp * tmp;
                }
                doubleTemp[i] = sum;
            }
        }
        
        else {
            for(i=0; i<numVectors; i++) {
                sum = 0.0;
                for(j=0; j<vnodeLength; j++) {
                    double tmp = values[i][j] / wValues[i][j];
                    sum += tmp * tmp;
                }
                doubleTemp[i] = sum;
            }
        }
        
        result = comm.sumAll(doubleTemp);
        double oneOverN = 1.0 / (double) globalLength;
        for(i=0; i<numVectors; i++) result[i] = Math.sqrt(result[i] * oneOverN);
        
        updateFlops(3 * globalLength * numVectors);
        
        return 0;
    }
    
    /**
     * Computes minimum value of each vector in calling multi-vector.
     *
     * @param result    result[i] will contain the minimum value of the ith vector
     */
    public int minValue(double [] result) {
        int i, j;
        for(i=0; i<numVectors; i++) {
            double minVal = values[i][0];
            for(j=1; j<vnodeLength; j++) minVal = Math.min(minVal, values[i][j]);
            doubleTemp[i] = minVal;
        }
        
        result = comm.minAll(doubleTemp);
        
        return 0;
    }
    
    /**
     * Computes maximum value of each vector in calling multi-vector.
     *
     * @param result    result[i] will contain the maximum value of the ith vector
     */
    public int maxValue(double [] result) {
        int i, j;
        for(i=0; i<numVectors; i++) {
            double maxVal = values[i][0];
            for(j=1; j<vnodeLength; j++) maxVal = Math.max(maxVal, values[i][j]);
            doubleTemp[i] = maxVal;
        }
        
        result = comm.maxAll(doubleTemp);
        
        return 0;
    }
    
    /**
     * Computes mean (average) value of each vector in multi-vector.
     *
     * @param result    result[i] contains the mean value of the ith vector
     */
    public int meanValue(double [] result) {
        int i, j;
        double fGlobalLength = 1.0 / Math.max((double)globalLength, 1.0);
        
        for(i=0; i<numVectors; i++) {
            double sum = 0.0;
            for(j=0; j<vnodeLength; j++) sum += values[i][j];
            doubleTemp[i] = sum;
        }
        
        result = comm.sumAll(doubleTemp);
        for(i=0; i<numVectors; i++) result[i] = result[i] * fGlobalLength;
        
        updateFlops(globalLength * numVectors);
        
        return 0;
    }
    
    /**
     * Performs matrix-matrix multiplication: this = scalar*this + scalarAB*a*b
     *
     * @param transA    "T" operates with transpose of a;
     *                  "N" operates with a (no transpose)
     *
     * @param transB    "T" operates with transpose of B;
     *                  "N" operates with b (no transpose)
     */
    public int multiply(String transA, String transB, double scalarAB,
    MultiVector a, MultiVector b, double scalar) {
        int aNRows = (transA.equals("T")) ? a.numVectors : a.vnodeLength;
        int aNCols = (transA.equals("T")) ? a.vnodeLength : a.numVectors;
        int bNRows = (transB.equals("T")) ? b.numVectors : b.vnodeLength;
        int bNCols = (transB.equals("T")) ? b.vnodeLength : b.numVectors;
        double scalarLocal = scalar;

        if(this.vnodeLength != aNRows ||
	   aNCols          != bNRows ||
	   this.numVectors != bNCols) return -2;
        
        boolean aIsLocal = !a.distributedGlobal;
        boolean bIsLocal = !b.distributedGlobal;
        boolean cIsLocal = !this.distributedGlobal;
        boolean case1 = ( aIsLocal &&  bIsLocal &&  cIsLocal);
        boolean case2 = (!aIsLocal && !bIsLocal &&  cIsLocal && transA == "T");
        boolean case3 = (!aIsLocal &&  bIsLocal && !cIsLocal && transA == "N");

        // Test for meaningful cases
        if(case1 || case2 || case3) {
            if(scalar != 0.0 && case2)
                if(comm.getVnodeID() != 0) scalarLocal = 0.0;

	    int m = this.vnodeLength;
	    int n = this.numVectors;
	    int k = aNCols;

	    double [][] aVals = a.values;
	    double [][] bVals = b.values;
	    double [][] cVals = this.values;

	    System.out.println("TransA: "+transA);
	    System.out.println("TransB: "+transB);
	    System.out.println("m:      "+m);
	    System.out.println("n:      "+n);
	    System.out.println("k:      "+k);
	    System.out.println("alpha:  "+scalarAB);
	    System.out.println("ColsA:  "+aVals.length);
	    System.out.println("RowsA:  "+aVals[0].length);
	    System.out.println("ColsB:  "+bVals.length);
	    System.out.println("RowsB:  "+bVals[0].length);
	    System.out.println("beta:   "+scalarLocal);
	    System.out.println("ColsC:  "+cVals.length);
	    System.out.println("RowsC:  "+cVals[0].length);

	    DGEMM.DGEMM(transA, transB, m, n, k, scalarAB, aVals, bVals,
			scalarLocal, cVals);
	    
	    if(case1) {
		updateFlops(2 * m * n * k);
		if(scalarAB != 1.0) updateFlops(m * n);
		if(scalar == 1.0) updateFlops(m * n);
		else if(scalar != 0.0) updateFlops(2 * m * n);
	    }
	    else if(case2) {
		updateFlops(2 * m * n * a.globalLength);
		if(scalarAB != 1.0) updateFlops(m * n);
		if(scalar == 1.0) updateFlops(m * n);
		else if(scalar != 0.0) updateFlops(2 * m * n);
	    }
	    else {
		updateFlops(2 * globalLength * n * k);
		if(scalarAB != 1.0) updateFlops(globalLength * n);
		if(scalar == 1.0) updateFlops(globalLength * n);
		else if(scalar != 0.0) updateFlops(2 * globalLength * n);
	    }
	    
	    if(case2) return(reduce());
	    return 0;

	}
	else return -3; // Not supported operation

    }
    
    /**
     * Multiplies a multi-vector with another element-by-element;
     * this = scalar*this + scalarAB*b@a where @ denotes element-wise multiplication
     */
    public int multiply(double scalarAB, MultiVector a, MultiVector b,
    double scalar) {
        int i, j;
        
        if(scalarAB == 0.0) return scale(scalar);
        if(a.numVectors != 1 && a.numVectors != b.numVectors) return -1;
        if(numVectors != b.numVectors) return -2;
        if(vnodeLength != a.vnodeLength || vnodeLength != b.vnodeLength) return -3;
        
        double [][] aVals = a.values;
        double [][] bVals = b.values;
        
        int incA = 1;
        if(a.numVectors == 1) incA = 0;
        if(scalar == 0.0) {
            if(scalarAB == 1.0) {
                for(i=0; i<numVectors; i++) {
                    double [] aVector = aVals[i*incA];
                    for(j=0; j<vnodeLength; j++) values[i][j] = aVector[j] * bVals[i][j];
                }
                updateFlops(globalLength * numVectors);
            }
                
            else {
                for(i=0; i<numVectors; i++) {
                    double [] aVector = aVals[i*incA];
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = scalarAB * aVector[j] * bVals[i][j];
                }
                updateFlops(2 * globalLength * numVectors);
            }
        }
        else if(scalar == 1.0) {
            if(scalarAB == 1.0) {
                for(i=0; i<numVectors; i++) {
                    double [] aVector = aVals[i*incA];
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] += aVector[j] * bVals[i][j];
                }
                updateFlops(2 * globalLength * numVectors);
            }
            else {
                for(i=0; i<numVectors; i++) {
                    double [] aVector = aVals[i*incA];
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] += scalarAB * aVector[j] * bVals[i][j];
                }
                updateFlops(3 * globalLength * numVectors);
            }
        }
        else { // scalar != 1.0 && scalar != 0
            if(scalarAB == 1.0) {
                for(i=0; i<numVectors; i++) {
                    double [] aVector = aVals[i*incA];
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = scalar * values[i][j] + aVector[j] * bVals[i][j];
                }
                updateFlops(3 * globalLength * numVectors);
            }
            else {
                for(i=0; i<numVectors; i++) {
                    double [] aVector = aVals[i*incA];
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = scalar * values[i][j] + scalarAB * aVector[j]
                        * bVals[i][j];
                }
                updateFlops(4 * globalLength * numVectors);
            }
        }
        return 0;
    }
    
    private int reduce() {
        // Global reduction on each entry of a replicated local multi-vector
        
        int i, j;
        
        double [] tmp = new double [vnodeLength];
        
        for(i=0; i<numVectors; i++) {
            for(j=0; j<vnodeLength; j++) tmp[j] = values[i][j];
            values[i] = comm.sumAll(tmp);
        }
        
        return 0;
}

    /**
     * Multiplies a multi-vector by the reciprical of another, element by element;
     * this = scalar*this + scalarAB*b@a where @ denotes element-wise division
     */
    public int reciprocalMultiply(double scalarAB, MultiVector a,
    MultiVector b, double scalar) {
        int i, j;
        
        if(scalarAB == 0.0) return scale(scalar);
        if(a.numVectors != 1 && a.numVectors != b.numVectors) return -1;
        if(numVectors != b.numVectors) return -2;
        if(vnodeLength != a.vnodeLength || vnodeLength != b.vnodeLength) return -3;
        
        double [][] aVals = a.values;
        double [][] bVals = b.values;
        
        int inCA = 1;
        if(a.numVectors == 1) inCA = 0;
        
        if(scalar == 0.0) {
            if(scalarAB == 1.0) {
                for(i=0; i<numVectors; i++) {
                    double [] aVector = aVals[i*inCA];
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = bVals[i][j] / aVector[j];
                }
                updateFlops(globalLength * numVectors);
            }
            else {
                for(i=0; i<numVectors; i++) {
                    double [] aVector = aVals[i*inCA];
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = scalarAB * bVals[i][j] / aVector[j];
                }
                updateFlops(2 * globalLength * numVectors);
            }
        }
        else if(scalar == 1.0) {
            if(scalarAB == 1.0) {
                for(i=0; i<numVectors; i++) {
                    double [] aVector = aVals[i*inCA];
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] += bVals[i][j] / aVector[j];
                }
                updateFlops(2 * globalLength * numVectors);
            }
            else {
                for(i=0; i<numVectors; i++) {
                    double [] aVector = aVals[i*inCA];
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] += scalarAB * bVals[i][j] / aVector[j];
                }
                updateFlops(3 * globalLength * numVectors);
            }
        }
        else { // scalar != 1.0 && scalar != 0
            if(scalarAB == 1.0) {
                for(i=0; i<numVectors; i++) {
                    double [] aVector = aVals[i*inCA];
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = scalar * values[i][j] + bVals[i][j] / aVector[j];
                }
                updateFlops(3 * globalLength * numVectors);
            }
            else {
                for(i=0; i<numVectors; i++) {
                    double [] aVector = aVals[i*inCA];
                    for(j=0; j<vnodeLength; j++)
                        values[i][j] = scalar * values[i][j] + scalarAB * bVals[i][j] / aVector[j];
                }
                updateFlops(4 * globalLength * numVectors);
            }
        }
        return 0;
    }
    
     /* DIFFERENT RANDOM FUNCTION
    public int setSeed(double seed) {
        this.seed = seed;
        return 0;
    }
    
    public double getSeed() {
        return seed;
    }
    */
    
    /**
     * Resets the number of floating point operations to zero.
     */
    public void resetFlops() {
        flops = 0.0;
    }
    
    /**
     * Accessor for numVectors.
     */
    public int getNumVectors() {
        return numVectors;
    }
    
    /**
     * Accessor for vnodeLength.
     */
    public int getVnodeLength() {
        return vnodeLength;
    }
    
    /**
     * Accessor for globalLength.
     */
    public int getGlobalLength() {
        return globalLength;
    }
    
    /**
     * Accessor for flops.
     */
    public double getFlops() {
        return flops;
    }
    
    /**
     * Accessor for map.
     */
    public BlockMap getMap() {
        return map;
    }
    
    /**
     * Accessor for comm.
     */
    public Comm getComm() {
        return comm;
    }
    
    /**
     * Accessor for distributeGlobal.
     */
    public boolean isDistributedGlobal() {
        return distributedGlobal;
    }
    
    protected double [][] getValues() {
        return values;
    }
    
    private void updateFlops(int flops) {
        this.flops = this.flops += (double) flops;
    }
    
    public double VALUEAT(int row, int col) {
        //if(row >= numVectors) throw new JpetraException("row must be less than numVectors");
        //if(col >= vnodeLength) throw new JpetraException("col must be less than vnodeLength");
        return values[row][col];
    }
    
    public Object clone() throws CloneNotSupportedException {
        MultiVector theClone = (MultiVector)super.clone();
        
        if(doubleTemp != null) {
            theClone.doubleTemp = new double [doubleTemp.length];
            System.arraycopy(doubleTemp, 0, theClone.doubleTemp, 0, doubleTemp.length);
        }
        else theClone.doubleTemp = null;
        
        if(exports != null) {
            theClone.exports = new double [exports.length];
            System.arraycopy(exports, 0, theClone.exports, 0, exports.length);
        }
        else theClone.exports = null;
        
        if(imports != null) {
            theClone.imports = new double [imports.length];
            System.arraycopy(imports, 0, theClone.imports, 0, imports.length);
        }
        else theClone.imports = null;
        
        if(intTemp != null) {
            theClone.intTemp = new int [intTemp.length];
            System.arraycopy(intTemp, 0, theClone.intTemp, 0, intTemp.length);
        }
        else theClone.intTemp = null;
        
        if(values != null) {
            theClone.values = new double [numVectors][vnodeLength];
            for(int i=0; i<numVectors; i++)
                System.arraycopy(values[i], 0, theClone.values[i], 0, values[i].length);
        }
        else theClone.values = null;
        
        return theClone;
    }
    public void print(String title) {
	System.out.println("\n\n*************************");
	System.out.println(title);
	System.out.println("*************************");
	print();
	System.out.println("*************************");
    }
    public void print() {
	System.out.println("Number Vectors:  "+getNumVectors());
	System.out.println("Global Length:   "+getGlobalLength());
	System.out.println("Process Length:     "+getVnodeLength());
	
	int numVecs = getNumVectors();
	int length = getVnodeLength();
	double [][] vals = this.extractView();
	
	System.out.println("   Process    Row Index    Col Index    Value");
	for(int i = 0; i<length; i++) {
	    int row = getMap().getGID(i);
	    for (int j=0; j<numVecs; j++)
		System.out.println("   1           "+row+"        "+j+"        "+vals[j][i]);
	}
    }
}

