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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

package Jpetra;

import java.util.Set;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import java.io.Serializable;
import java.io.Externalizable;
import java.io.ObjectOutput;
import java.io.ObjectInput;
import java.io.ObjectOutputStream;
import java.io.ObjectInputStream;
import java.io.FileOutputStream;
import java.io.FileInputStream;

/**
 *
 * @author  Jason Cross
 */
public class CisMatrix extends DistObject implements Externalizable {
    /**
     * Pass to the CisMatrix constructor to make the CisMatrix row oriented.
     */
    public static final boolean ROW_ORIENTED = true;
    /**
     * Pass to the CisMatrix constructor to make the CisMatrix col oriented.
     */
    public static final boolean COL_ORIENTED = false;
    
    public static final boolean USE_TRANSPOSE_A = true;
    public static final boolean USE_A = false;
    
    private Graph graph;  // nonzero pattern of the matrix;  secondary indices for doubleValues
    private VectorSpace primaryVectorSpace;
    private VectorSpace secondaryVectorSpace;
    private boolean rowOriented;  // == either ROW_ORIENTED or COL_ORIENTED
    private boolean filled;  // true if fillComplete() has been called an no methods modifying the matrix have been called
    private double[] doubleValues;  // the actual values of the matrix
    private int[] numEntries;  // number of entries per row/col
    private int[] startIndex;   // number of entries before the current row/col
    private TreeMap outerTree;  // map of maps used to construct the matrix in its unoptimized form
    private TreeSet secondaryTree; // map of nonzero vectors which is used to construct the secondaryVectorSpace
    private int numTotalEntries;  // total number of values in the matrix
    /*private int maxSecondaryId;  // used to build the secondaryVectorSpace automatically, may not be needed!!!*/
    
    // for reverse op
    private boolean doneForward;  // true if the forward import/export have been performed on this CisMatrix
    
    public CisMatrix() {
        this.secondaryTree = new TreeSet();
    }
    
    /**
     * Construct an empty <code>CisMatrix</code>.  The <code>secondaryVectorSpace</code> will be generated automatically
     * once <code>fillComplete</code> is called.
     *
     * @param primaryVectorSpace describes what global elemenets this <code>CisMatrix</code> owns
     * @param rowOriented determines whether the <code>CisMatrix</code> is row or column oriented
     */
    public CisMatrix(VectorSpace primaryVectorSpace, boolean rowOriented) {
        //this.outputStreams.put("CISMATRIX", new Output("CisMatrix: ", true, System.out, false, System.out));
        
        this.filled = false;
        /*this.maxSecondaryId = 0;*/
        this.primaryVectorSpace=primaryVectorSpace;
        this.rowOriented = rowOriented;
        this.outerTree = new TreeMap();
        this.secondaryTree = new TreeSet();
    }
    
    public void insertEntries(int localRowColId, int[] indices, double[] entries, int combineMode) {
        /*if (this.filled) {
            this.println("ERR", "insertEntiries cannot be called after fillComplete().");
        }*/
        this.filled = false;
        
        // !! need to check if I own the globalRowColId
        
        // need to see if the TreeMap exists for the specified row/col
        TreeMap rowColTreeMap;
        if (!this.outerTree.containsKey(new Integer(localRowColId))) {
            rowColTreeMap = new TreeMap();
            //this.println("STD", "globalRowCol does not exist, creating...");
            this.outerTree.put(new Integer(localRowColId), rowColTreeMap);
        }
        else {
            //this.println("STD", "localRowCol exists, setting rowColTreeMap to existing TreeMap...");
            rowColTreeMap = (TreeMap) outerTree.get(new Integer(localRowColId));
        }
        
        // now that we know the row/col exists, insert entries into the row/col
        // and sum them with any pre-existing entries
        //this.println("STD", "CisMatrix is internally inserting entries...");
        Object temp;
        int toInsert;
        if (combineMode == DistObject.ADD) {
            for(int i=0; i < indices.length; i++) {
                //this.println("STD", "Inserting index " + indices[i]);
                temp = rowColTreeMap.get(new Integer(indices[i]));
                if (temp == null) {
                    //this.println("STD", "Index " + indices[i] + " does not exist, creating...");
                    this.numTotalEntries++;
                    rowColTreeMap.put(new Integer(indices[i]), new Double(entries[i]));
                }
                else {
                    //this.println("STD", "Index " + indices[i] + " does exists, adding to prevous value...");
                    rowColTreeMap.put(new Integer(indices[i]), new Double(entries[i] + ((Double) temp).doubleValue()));
                }
            }
        }
        else if (combineMode == DistObject.REPLACE) {
            for(int i=0; i < indices.length; i++) {
                this.numTotalEntries++;
                rowColTreeMap.put(new Integer(indices[i]), new Double(entries[i]));
            }
        }
        else {
            this.println("ERR", "The specified combine mode is not supported by CisMatrix insertEntries.");
            System.exit(1);
        }
    }
    
    public void insertEntry(int localRowColId, int index, double entry, int combineMode) {
        /*if (this.filled) {
            this.println("ERR", "insertEntiry cannot be called after fillComplete().");
        }*/
        this.filled = false;
        
        // !! need to check if I own the globalRowColId
        
        // need to see if the TreeMap exists for the specified row/col
        TreeMap rowColTreeMap;
        if (!this.outerTree.containsKey(new Integer(localRowColId))) {
            rowColTreeMap = new TreeMap();
            //this.println("STD", "globalRowCol does not exist, creating...");
            this.outerTree.put(new Integer(localRowColId), rowColTreeMap);
        }
        else {
            //this.println("STD", "globalRowCol exists, setting rowColTreeMap to existing TreeMap...");
            rowColTreeMap = (TreeMap) outerTree.get(new Integer(localRowColId));
        }
        
        // now that we know the row/col exists, insert entries into the row/col
        // and sum them with any pre-existing entries
        //this.println("STD", "CisMatrix is internally inserting an entry...");
        Object temp;
        int toInsert;
        //this.println("STD", "Inserting index " + index);
        if (combineMode == DistObject.ADD) {
            temp = rowColTreeMap.get(new Integer(index));
            if (temp == null) {
                //this.println("STD", "Index " + index + " does not exist, creating...");
                this.numTotalEntries++;
                rowColTreeMap.put(new Integer(index), new Double(entry));
            }
            else {
                //this.println("STD", "Index " + index + " does exists, adding to prevous value...");
                rowColTreeMap.put(new Integer(index), new Double(entry + ((Double) temp).doubleValue()));
            }
        } else if (combineMode == DistObject.REPLACE) {
            this.numTotalEntries++;
            rowColTreeMap.put(new Integer(index), new Double(entry));
        }
        else {
            this.println("ERR", "The specified combine mode is not supported by CisMatrix insertEntry.");
            System.exit(1);
        }
    }
    
    public void fillComplete() {
        /*if (this.filled) {
            this.print("ERR", "fillComplete() may only be called once.");
        }*/
        this.filled = true;
        this.doubleValues = new double[this.numTotalEntries];
        this.numEntries = new int[primaryVectorSpace.getNumMyGlobalEntries()];
        this.startIndex = new int[primaryVectorSpace.getNumMyGlobalEntries()];
        
        Set outerKeysValues = this.outerTree.entrySet();
        Iterator outerIterator = outerKeysValues.iterator();
        Map.Entry outerMapEntry;
        
        TreeMap innerTree;
        Set innerKeysValues;
        Iterator innerIterator;
        Map.Entry innerMapEntry;
        
        int nextEntryIndex=0;
        int numEntriesColRow;
        int[] tempGraph = new int[this.numTotalEntries];
        int i=0;
        int intStartIndex = 0;
        Integer secondaryIntegerId;
        while(outerIterator.hasNext()) {
            //this.println("STD", "Doing an outer loop...");
            outerMapEntry = (Map.Entry) outerIterator.next();
            innerTree = (TreeMap) outerMapEntry.getValue();
            innerKeysValues = innerTree.entrySet();
            innerIterator = innerKeysValues.iterator();
            numEntriesColRow = 0;
            while (innerIterator.hasNext()) {
                //this.println("STD", "Doing an inner loop...");
                innerMapEntry = (Map.Entry) innerIterator.next();
                this.doubleValues[nextEntryIndex] = ((Double) innerMapEntry.getValue()).doubleValue();
                secondaryIntegerId = (Integer) innerMapEntry.getKey();
                if (!this.secondaryTree.contains(secondaryIntegerId)) {
                    this.secondaryTree.add(secondaryIntegerId);
                }
                tempGraph[nextEntryIndex++] = secondaryIntegerId.intValue();
                numEntriesColRow++;
            }
            /*if (this.maxSecondaryId < tempGraph[nextEntryIndex-1]) {
                this.maxSecondaryId = tempGraph[nextEntryIndex-1];
            }*/
            this.startIndex[i] = intStartIndex;
            intStartIndex += numEntriesColRow;
            this.numEntries[i++] = numEntriesColRow;
        }
        this.graph = new Graph(primaryVectorSpace, tempGraph, this.numEntries);
        
        // now build the secondaryVectorSpace from secondaryTree
        Iterator secondaryIterator = secondaryTree.iterator();
        int[] secondaryGids = new int[secondaryTree.size()];
        i=0;
        int tmp;
        while (secondaryIterator.hasNext()) {
            secondaryGids[i++] = ((Integer) secondaryIterator.next()).intValue();
        }
        this.secondaryVectorSpace = new VectorSpace(new ElementSpace(secondaryGids, this.primaryVectorSpace.getComm()));
        
        // if we wanted to do away with the TreeMaps at this point to save memory
        // then outerTree map would need to be set to null so the java garbage collector
        // would know to pick it up
        // outerTree = null;
        // secondaryTree = null;
    }
    
    public void scale(double scaler) {
        if (!this.filled) {
            this.println("FATALERR", "You must call fillComplete() before calling scale(double scaler).");
            System.exit(0);
        }
        
        for(int i=0; i < this.doubleValues.length; i++) {
            this.doubleValues[i] *= scaler;
        }
    }
    
    public void printOutAllVnodes(String iostream) {
        this.printOut(iostream, true);
    }
    
    public void printOut(String iostream) {
        this.printOut(iostream, false);
    }
    
    public void printOut(String iostream, boolean all) {
        if (all || (primaryVectorSpace.getComm().getVnodeId() == 0)) {
            this.println("STD", "CisMatrix.printOut() is starting...");
            if (!filled) {
                this.print("ERR", "You must call fillComplete() before calling printOut(String iostream)");
                return;
            }
            
            int i=0;
            if (this.rowOriented) {
                //this.println("STD", "this.numEntries.length=" + numEntries.length);
                for(int row=0; row < this.numEntries.length; row++) {
                    //this.println("STD", "Doing a row loop...");
                    for(int entry=0; entry < numEntries[row]; entry++) {
                        //this.println("STD", "\nDoing an entry loop...");
                        this.println(iostream, row + " " + graph.getIndex(i) + " " + doubleValues[i++]);
                    }
                }
            }
            else {
                for(int col=0; col < this.numEntries.length; col++) {
                    //this.println("STD", "Doing a row loop...");
                    for(int entry=0; entry < numEntries[col]; entry++) {
                        //this.println("STD", "\nDoing an entry loop...");
                        this.println(iostream, graph.getIndex(i) + " " + col + " " + doubleValues[i++]);
                    }
                }
                
            }
            this.println("STD", "CisMatrix.printOut() has ended...");
        }
    }
    
    public int getNumMyNonZeros() {
        return this.doubleValues.length;
    }
    
    public int getNumGlobalRows() {
        if (this.rowOriented) {
            return this.primaryVectorSpace.getNumGlobalEntries();
        }
        else {
            return this.secondaryVectorSpace.getNumGlobalEntries();
        }
    }
    
    public int getNumGlobalColumns() {
        if (this.rowOriented) {
            return this.secondaryVectorSpace.getNumGlobalEntries();
        }
        else {
            return this.primaryVectorSpace.getNumGlobalEntries();
        }
    }
    
    public int getNumMyRows() {
        if (this.rowOriented) {
            return this.primaryVectorSpace.getNumMyGlobalEntries();
        }
        else {
            return this.secondaryVectorSpace.getNumMyGlobalEntries();
        }
    }
    
    public int getNumMyColumns() {
        if (this.rowOriented) {
            return this.secondaryVectorSpace.getNumMyGlobalEntries();
        }
        else {
            return this.primaryVectorSpace.getNumMyGlobalEntries();
        }
    }
    
    public boolean isRowOriented() {
        return this.rowOriented;
    }
    
    public Graph getGraph() {
        return this.graph;
    }
    
    public double[] getEntriesArray() {
        return this.doubleValues;
    }
    
    public int[] getNumEntriesArray() {
        return this.numEntries;
    }
    
    public int[] getStartIndex() {
        return this.startIndex;
    }
    
    public int[] getNumEntries() {
        return this.numEntries;
    }
    
    public boolean isFilled() {
        return this.filled;
    }
    
    public VectorSpace getVectorSpace() {
        return this.primaryVectorSpace;
    }
    
    public Serializable[] packAndPrepare(DistObject distObjectSource, int[] exportGids, int[] exportLids) {
        CisMatrix exportMatrix = (CisMatrix) distObjectSource;
        
        Serializable[] exportData = new Serializable[exportGids.length];  // the object to be exported by Distributor.distribute
        
        // this code depends on this CisMatrix being filled
        if (!exportMatrix.isFilled()) {
            exportMatrix.fillComplete();
        }
        
        // get all necessary information to allow us to access the elements of exportMatrix
        int[] exportNumEntries = exportMatrix.getNumEntries();
        int[] exportStartIndex = exportMatrix.getStartIndex();
        double[] exportDoubleValues = exportMatrix.getEntriesArray();
        int[] graphArray = exportMatrix.getGraph().getNonZeroEntriesArray();
        
        int[] indices;
        double[] values;
        Serializable[] element;
        int lid;
        Serializable[] noEntries = new Serializable[]{new Integer(-1)};  // send this array if the entire elements (row/col) of exportMatrix is filled with 0's
        for(int i=0; i < exportLids.length; i++) {
            lid = exportLids[i];
            if (exportNumEntries[lid] == 0) {
                exportData[i] = noEntries;
            }
            else {
                // copy array slices from the exportDoubleValues and graphArray into temporary arrays
                // so that the element (col/row) can be sent to another vnode
                indices = new int[exportNumEntries[lid]];
                values = new double[exportNumEntries[lid]];
                System.arraycopy(exportDoubleValues, exportStartIndex[lid], values, 0, values.length);
                System.arraycopy(graphArray, exportStartIndex[lid], indices, 0, indices.length);
                
                // data is prepared, now pack it up
                element = new Serializable[3];
                element[0] = new Integer(exportGids[i]);
                element[1] = indices;
                element[2] = values;
                exportData[i] = element;
            }
        }
        
        return exportData;
    }
    
    public int[][] unpackAndCombine(Serializable[] importData, int combineMode) {
        int[] indices;
        double[] values;
        Serializable[] element;
        Serializable[] elementArray;
        int gid;
        int lid;
        
        // for reverse op
        int[][] reverseExportVnodeIdsGidsLids = null;
        if (!this.doneForward) {
            int sumEntries = 0;
            for(int i=0; i < importData.length; i++) {
                // if a vnode didn't send us any data, then importData[vnodeId] == null
                if (importData[i] != null) {
                    sumEntries += ((Serializable[]) importData[i]).length;
                }
            }
            reverseExportVnodeIdsGidsLids = new int[4][];
            reverseExportVnodeIdsGidsLids[0] = new int[primaryVectorSpace.getComm().getNumVnodes()];
            reverseExportVnodeIdsGidsLids[1] = new int[sumEntries];
            reverseExportVnodeIdsGidsLids[2] = new int[sumEntries];
            reverseExportVnodeIdsGidsLids[3] = new int[sumEntries];
        }
        int revCount = 0;
        for(int i=0; i < importData.length; i++) {
            if (importData[i] != null) {
                elementArray = (Serializable[]) importData[i];
                for(int j=0; j < elementArray.length; j++) {
                    element = (Serializable[]) elementArray[j];
                    gid = ((Integer) element[0]).intValue();
                    //this.println("CISMATRIX", "Checking gid " + gid + " from vnode " + i);
                    lid = primaryVectorSpace.getLocalIndex(gid);
                    // gid == -1 means that the sending vnode didn't have a nonzero value for the entire
                    // row/col of that gid so we can just ignore it
                    if (!this.doneForward) {
                        reverseExportVnodeIdsGidsLids[0][i]++;
                        reverseExportVnodeIdsGidsLids[1][revCount] = i;
                        reverseExportVnodeIdsGidsLids[2][revCount] = gid;
                        reverseExportVnodeIdsGidsLids[3][revCount++] = lid;
                    }
                    if (gid != -1 && lid != -1) {
                        //this.println("CISMATRIX", "adding " + ((int[]) element[1]).length + " elements to gid " + gid + " from vnode " + i);
                        this.insertEntries(lid, (int[]) element[1], (double[]) element[2], combineMode);
                    }
                } // end for j
            } // end if (importData[i] != null)
        }  // end for i
        
        if (!this.doneForward) {
            this.doneForward = true;
        }
        
        return reverseExportVnodeIdsGidsLids;
    }
    
    public void copyAndPermute(DistObject distObjectSource, int numSameGids, int[] permuteToLids, int[] permuteFromLids, int combineMode) {
        CisMatrix sourceMatrix = (CisMatrix) distObjectSource;
        // get all the necessary information from sourceMatrix in order to access its elements
        double[] srcValues = sourceMatrix.getEntriesArray();
        int[] srcStartIndices = sourceMatrix.getStartIndex();
        int[] srcIndices = (sourceMatrix.getGraph()).getNonZeroEntriesArray();
        int[] srcNumEntries = sourceMatrix.getNumEntries();
        
        double[] tmpEntries;
        int[] tmpIndices;
        for(int i=0; i < numSameGids; i++) {
            // srcNumEntries[i] == 0 then all the values in that row/col are 0 so we don't need to do anything
            if (srcNumEntries[i] != 0) {
                // copy array slices from the srcValues and srcIndices into temporary arrays
                // so that insertEntries can be used to insert the imported data into this CisMatrix
                tmpEntries = new double[srcNumEntries[i]];
                tmpIndices = new int[srcNumEntries[i]];
                System.arraycopy(srcValues, srcStartIndices[i], tmpEntries, 0, tmpEntries.length);
                System.arraycopy(srcIndices, srcStartIndices[i], tmpIndices, 0, tmpIndices.length);
                this.insertEntries(i, tmpIndices, tmpEntries, combineMode);
            }
        }
        
        int srcLid;
        for(int i=0; i < permuteToLids.length; i++) {
            srcLid = permuteFromLids[i];
            if (srcNumEntries[srcLid] != 0) {
                // this does the same thing as the for loop above with numSameGids except
                // here we also permute the values from the src Gid to the target Gid
                tmpEntries = new double[srcNumEntries[srcLid]];
                tmpIndices = new int[srcNumEntries[srcLid]];
                System.arraycopy(srcValues, srcStartIndices[srcLid], tmpEntries, 0, tmpEntries.length);
                System.arraycopy(srcIndices, srcStartIndices[srcLid], tmpIndices, 0, tmpIndices.length);
                this.insertEntries(permuteToLids[i], tmpIndices, tmpEntries, combineMode);
            }
        }
    }
    
    public VectorSpace getColumnVectorSpace() {
        if (this.rowOriented) {
            if (this.filled) {
                return this.secondaryVectorSpace;
            }
            else {
                this.println("FATALERR", "You cannot retrive the secondaryVectorSpace (in this case the columnVectorSpace) before you call fillComplete().");
                System.exit(1);
            }
        }
        
        return this.primaryVectorSpace;
    }
    
    public VectorSpace getRowVectorSpace() {
        if (this.rowOriented) {
            return this.primaryVectorSpace;
        }
        
        if (this.filled) {
            return this.secondaryVectorSpace;
        }
        else {
            this.println("FATALERR", "You cannot retrive the secondaryVectorSpace (in this case the rowVectorSpace) before you call fillComplete().");
            System.exit(1);
        }
        
        return null;  // will never get here, but java complains if this is not present
    }
    
    public VectorSpace getPrimaryVectorSpace() {
        return this.primaryVectorSpace;
    }
    
    public VectorSpace getSecondaryVectorSpace() {
        if (this.filled) {
            return this.secondaryVectorSpace;
        }
        else {
            this.println("FATALERR", "You cannot retrive the secondaryVectorSpace before you call fillComplete().");
            System.exit(1);
        }
        
        return null;  // will never get here, but java complains if this is not present
    }
    
    public void multiply(boolean useTransposeA, MultiVector x, MultiVector y) {
        Time timer = new Time(this.primaryVectorSpace.getComm());
        if (!this.filled) {
            this.println("FATALERR", "A CisMatrix must be filled before multiply can be called on it.");
            System.exit(1);
        }
        
        // zero out y
        y.putScalar(0.0);
        
        // setup temporary MultiVectors
        
        /*
        int[] tmp = x.getVectorSpace().getMyGlobalEntryIds();
        for(int i=0; i < tmp.length; i++) {
            this.println("CISMATRIX", "my x Gid: " + tmp[i]);
        }
        tmp = importMultiVector.getVectorSpace().getMyGlobalEntryIds();
        for(int i=0; i < tmp.length; i++) {
            this.println("CISMATRIX", "my importMultiVector Gid: " + tmp[i]);
        }
         */
        
        // import any values needed for x from other vnodes
        
        // compatibility checks
        //if (useTransposeA == CisMatrix.USE_A) {
        if (useTransposeA == CisMatrix.USE_A) {
            //this.println("CISMATRIX", "Doing compatibility checks for A.");
            if (this.getColumnVectorSpace().getNumGlobalEntries() > x.getVectorSpace().getNumGlobalEntries()) {
                this.println("FATALERR", "In CisMatrix.muliptly: The number of columns in CisMatrix A (" + this.getColumnVectorSpace().getNumGlobalEntries() + ") > the number of rows (" + x.getVectorSpace().getNumGlobalEntries() + ") in MultiVector x.");
                System.exit(1);
            }
            if (this.getRowVectorSpace().getNumGlobalEntries() > y.getVectorSpace().getNumGlobalEntries()) {
                this.println("FATALERR", "In CisMatrix.muliptly: The number of rows in CisMatrix A (" + this.getRowVectorSpace().getNumGlobalEntries() + ") > the number of rows (" + y.getVectorSpace().getNumGlobalEntries() + ") in MultiVector y.");
                System.exit(1);
            }
            if (x.getNumCols() != y.getNumCols()) {
                this.println("FATALERR", "In CisMatrix.muliptly: The number of columns in MultiVector x (" + x.getNumCols() + ") != the number of columns (" + y.getNumCols() + ") in MultiVector y.");
                System.exit(1);
            }
        }
        else {
            this.println("CISMATRIX", "Doing compatibility checks for the transpose of A.");
            if (this.getRowVectorSpace().getNumGlobalEntries() > x.getVectorSpace().getNumGlobalEntries()) {
                this.println("FATALERR", "In CisMatrix.multiply: The number of columns in CisMatrix A' (" + this.getRowVectorSpace().getNumGlobalEntries() + ") > the number of rows (" + x.getVectorSpace().getNumGlobalEntries() + ") in MultiVector x.");
                System.exit(1);
            }
            if (this.getColumnVectorSpace().getNumGlobalEntries() > y.getVectorSpace().getNumGlobalEntries()) {
                this.println("FATALERR", "In CisMatrix.multiply: The number of rows in CisMatrix A' (" + this.getColumnVectorSpace().getNumGlobalEntries() + ") > the number of rows (" + y.getVectorSpace().getNumGlobalEntries() + ") in MultiVector y.");
                System.exit(1);
            }
            if (x.getNumCols() != y.getNumCols()) {
                this.println("FATALERR", "In CisMatrix.multiply: The number of columns in MultiVector x (" + x.getNumCols() + ") != the number of columns (" + y.getNumCols() + ") in MultiVector y.");
                System.exit(1);
            }
        }
        
        double[][] exportValues;
        MultiVector exportMultiVector;
        
        int[] indices = this.graph.getNonZeroEntriesArray();
        double sum;
        int index;
        if(useTransposeA == CisMatrix.USE_A) {
            if(this.rowOriented) {
                MultiVector importMultiVector = new MultiVector(this.getColumnVectorSpace(), new double[x.getNumCols()][this.getNumGlobalColumns()]);
                Import importer = new Import(x.getVectorSpace(), importMultiVector.getVectorSpace());
                importMultiVector.importValues(x, importer, DistObject.REPLACE);
                //this.println("CISMATRIX", "printing importMultiVector...");
                //importMultiVector.printOutAllVnodes("STD");
                double[][] importValues = importMultiVector.getValues();
                
                exportValues = new double[importValues.length][this.getNumMyRows()];
                exportMultiVector = new MultiVector(this.getRowVectorSpace(), exportValues);
                VectorSpace importMultiVectorVectorSpace = importMultiVector.getVectorSpace();
                double setupTime = timer.getElapsedTime();
                timer.resetStartTime();
                for(int row=0; row < this.getNumMyRows(); row++) {
                    for(int vector=0; vector < importMultiVector.getNumCols(); vector++){
                        sum = 0;
                        index = this.startIndex[row];
                        for(int col=0; col < this.numEntries[row]; col++) {
                            //this.println("STD", "this.doubleValues[index] * importValues[vector][indices[index++]]" + this.doubleValues[index] +" * " + importValues[vector][importMultiVectorVectorSpace.getLocalIndex(indices[index])]); //importValues[vector][indices[index]]);
                            //sum += this.doubleValues[index] * importValues[vector][importMultiVectorVectorSpace.getLocalIndex(indices[index++])];
			    
			    // NOTE!! i think importValues[vector][col] is incorrect and should be importValues[vector][importMultiVectorVectorSpace.getLocalIndex(indices[index++])]
                            //sum += this.doubleValues[index++] * importValues[vector][col];  // assume that Gids of importValues are in the same order as the Gids from doubleValues
                            sum += this.doubleValues[index] * importValues[vector][importMultiVectorVectorSpace.getLocalIndex(indices[index++])];
                        }
                        exportValues[vector][row] = sum;
                    }
                }
                double multTime = timer.getElapsedTime();
                timer.resetStartTime();
                Export exporter = new Export(exportMultiVector.getVectorSpace(), y.getVectorSpace());
                y.exportValues(exportMultiVector, exporter, DistObject.ADD);
                double exportTime = timer.getElapsedTime();
                if (this.primaryVectorSpace.getComm().getNumVnodes() == 1) {
                    this.println("STD", "setupTime\t" + setupTime + "\nmultTime\t" + multTime + "\nexportTime\t" + exportTime);
                }
                else {
                    this.println("STD", setupTime + "\n" + multTime + "\n" + exportTime);
                }
            }
            else {
                // possibly generate a new VectorSpace for the importMultiVector from the nonzero vectors of the ColumnVectorSpace
                // would be less communication if there are 0 vectors
                // for now this works though
                MultiVector importMultiVector = new MultiVector(this.getColumnVectorSpace(), new double[x.getNumCols()][this.getNumMyColumns()]);
                Import importer = new Import(x.getVectorSpace(), importMultiVector.getVectorSpace());
                importMultiVector.importValues(x, importer, DistObject.REPLACE);
                double[][] importValues = importMultiVector.getValues();
                exportValues = new double[importValues.length][this.getNumMyRows()];
                exportMultiVector = new MultiVector(this.getRowVectorSpace(), exportValues);
                //this.println("CISMATRIX", "y.getVectorSpace().getNumGlobalEntries(): " + y.getVectorSpace().getNumGlobalEntries());
                //this.println("CISMATRIX", "exportMultiVector.getVectorSpace().getNumGlobalEntries(): " + exportMultiVector.getVectorSpace().getNumGlobalEntries());
                for(int col=0; col < this.getNumMyColumns(); col++) {
                    for(int vector=0; vector < importMultiVector.getNumCols(); vector++){
                        index = this.startIndex[col];
                        for(int row=0; row < this.numEntries[col]; row++) {
                            exportValues[vector][row] += this.doubleValues[index++] * importValues[vector][col];
                        }
                    }
                }
                Import import2 = new Import(y.getVectorSpace(), exportMultiVector.getVectorSpace());
                y.exportValues(exportMultiVector, import2, DistObject.ADD);
            }
        }
        else {
            if(this.rowOriented) {
                exportValues = new double[x.getNumCols()][this.getNumMyColumns()];
                exportMultiVector = new MultiVector(this.getColumnVectorSpace(), exportValues);
                Import importer = new Import(x.getVectorSpace(), this.getColumnVectorSpace());
                exportMultiVector.importValues(x, importer, DistObject.REPLACE);
                // need to check this....
                MultiVector importMultiVector = new MultiVector(this.getColumnVectorSpace(), new double[x.getNumCols()][this.getNumMyColumns()]);
                VectorSpace importMultiVectorVS = importMultiVector.getVectorSpace();
                double[][] importValues = importMultiVector.getValues();
                this.fillComplete();
                //this.println("CISMATRIX", "doubleValues.length: " + doubleValues.length);
                //this.printOutAllVnodes("CISMATRIX");
                
                for(int row=0; row < this.getNumMyRows(); row++) {
                    for(int vector=0; vector < exportMultiVector.getNumCols(); vector++){
                        index = this.startIndex[row];
                        for(int col=0; col < this.numEntries[row]; col++) {
                            //importValues[vector][indices[index]] += this.doubleValues[index++] * exportValues[vector][row];
                            importValues[vector][importMultiVectorVS.getLocalIndex(indices[index])] += this.doubleValues[index++] * exportValues[vector][row];
                        }
                    }
                }
                
                Import import2 = new Import(y.getVectorSpace(), importMultiVector.getVectorSpace());
                y.exportValues(importMultiVector, import2, DistObject.ADD);
                //Export exporter = new Export(y.getVectorSpace(), importMultiVector.getVectorSpace());
                //y.exportValues(importMultiVector, exporter, DistObject.ADD);
            }
            else {
                MultiVector importMultiVector = new MultiVector(this.getRowVectorSpace(), new double[x.getNumCols()][this.getNumMyRows()]);
                Import importer = new Import(x.getVectorSpace(), importMultiVector.getVectorSpace());
                importMultiVector.importValues(x, importer, DistObject.REPLACE);
                double[][] importValues = importMultiVector.getValues();
                
                exportValues = new double[importValues.length][this.getNumMyColumns()];
                exportMultiVector = new MultiVector(this.getColumnVectorSpace(), exportValues);
                for(int col=0; col < this.getNumMyColumns(); col++) {
                    for(int vector=0; vector < importMultiVector.getNumCols(); vector++){
                        sum = 0;
                        index = this.startIndex[col];
                        for(int row=0; row < this.numEntries[col]; row++) {
                            sum += this.doubleValues[index] * importValues[vector][indices[index++]];
                        }
                        exportValues[vector][col] = sum;
                    }
                }
                Export exporter = new Export(exportMultiVector.getVectorSpace(), y.getVectorSpace());
                y.exportValues(exportMultiVector, exporter, DistObject.ADD);
                //Import import2 = new Import(y.getVectorSpace(), exportMultiVector.getVectorSpace());
                //y.exportValues(exportMultiVector, import2, DistObject.ADD);
            }
        }
        
        this.updateFlops(x.getNumCols() * this.getNumMyNonZeros() * 2.0);
    }
    
    public void writeExternal(ObjectOutput out) throws java.io.IOException {
        out.writeObject(this.primaryVectorSpace);
        out.writeBoolean(this.rowOriented);
        out.writeObject(this.doubleValues);
        out.writeObject(this.numEntries);
        out.writeObject(this.startIndex);
        out.writeObject(this.outerTree);
        out.writeInt(this.numTotalEntries);
    }
    
    public void readExternal(ObjectInput in) throws java.io.IOException, ClassNotFoundException {
        this.primaryVectorSpace = (VectorSpace) in.readObject();
        this.rowOriented = in.readBoolean();
        this.doubleValues = (double[]) in.readObject();
        this.numEntries = (int[]) in.readObject();
        this.startIndex = (int[]) in.readObject();
        this.outerTree = (TreeMap) in.readObject();
        this.numTotalEntries = in.readInt();
    }
    
    public static CisMatrix readFromFile(String fileName, Comm comm) {
        if (comm.getVnodeId() == 0) {
            try {
                ObjectInputStream ois = new ObjectInputStream(new FileInputStream(fileName));
                CisMatrix cisMatrix = (CisMatrix) ois.readObject();
                ois.close();
                cisMatrix.getVectorSpace().setComm(comm);
                return cisMatrix;
            } catch (java.io.IOException e) {
                JpetraObject.println("FATALERR", e.toString());
                System.exit(1);
            } catch (java.lang.ClassNotFoundException e) {
                JpetraObject.println("FATALERR", e.toString());
                System.exit(1);
            }
        }
        
        return null;
    }
    
    public static void writeToFile(String fileName, CisMatrix cisMatrix) {
        if (cisMatrix.primaryVectorSpace.getComm().getVnodeId() == 0) {
            if (cisMatrix.getVectorSpace().isDistributedGlobally()) {
                JpetraObject.println("FATALERR", "Only a serial MultiVector can be written to a serialized object file.");
                System.exit(1);
            }
            
            try {
                ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(fileName));
                oos.writeObject(cisMatrix);
                oos.close();
            } catch (java.io.IOException e) {
                JpetraObject.println("FATALERR", e.toString());
                System.exit(1);
            }
        }
    }
    
    public Object clone() {
        CisMatrix cloneCisMatrix = new CisMatrix(this.primaryVectorSpace, this.rowOriented);
        cloneCisMatrix.outerTree = (TreeMap) this.outerTree.clone();
        cloneCisMatrix.numTotalEntries = this.numTotalEntries;
        cloneCisMatrix.setFlopCounter(this.getFlopCounter());
        
        return cloneCisMatrix;
    }
    
    public boolean equals(Object obj) {
        // check to see if the reference addresses are the same
        if (this == obj) {
            return true;
        }
        
        CisMatrix otherCisMatrix = (CisMatrix) obj;
        if (!otherCisMatrix.primaryVectorSpace.equals(this.primaryVectorSpace)) {
            return false;
        }
        
        if (!otherCisMatrix.outerTree.equals(this.outerTree)) {
            return false;
        }
        
        return true;
    }
}
