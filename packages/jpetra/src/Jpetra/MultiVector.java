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

package Jpetra;

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
public class MultiVector extends DistObject implements Externalizable {
    VectorSpace vectorSpace;
    double[][] values;  // in column major form
    private boolean doneForward;
    
    public MultiVector() {
        // empty
    }
    
    public MultiVector(VectorSpace vectorSpace) {
        this.vectorSpace = vectorSpace;
    }
    
    public MultiVector(VectorSpace vectorSpace, double[][] values) {
        this.vectorSpace = vectorSpace;
        this.values = values;
    }
    
    public void putScalar(double scalar) {
        for(int i=0; i < this.values.length; i++) {
            for (int j=0; j < this.values[i].length; j++) {
                this.values[i][j] = scalar;
            }
        }
    }
    
    public double[] maxValue() {
        double maxValue;
        double[] result = new double[values.length];
        for(int i=0; i < values.length; i++) {
            maxValue = values[i][0];
            for(int j=1; j < values[i].length; j++) {
                maxValue = Util.max(maxValue, values[i][j]);
            }
            result[i] = maxValue;
        }
        
        result = vectorSpace.getComm().maxAll(result);
        return result;
    }
    
    public double[] minValue() {
        double minValue;
        double[] result = new double[values.length];
        for(int i=0; i < values.length; i++) {
            minValue = values[i][0];
            for(int j=1; j < values[i].length; j++) {
                minValue = Util.min(minValue, values[i][j]);
            }
            result[i] = minValue;
        }
        
        result = vectorSpace.getComm().minAll(result);
        return result;
    }
    
    public double[] meanValue() {
        double sum;
        double[] result = new double[values.length];
        for(int i=0; i < values.length; i++) {
            sum = values[i][0];
            for(int j=1; j < values[i].length; j++) {
                sum += values[i][j];
            }
            result[i] = sum;
        }
        
        result = vectorSpace.getComm().sumAll(result);
        int numGlobalEntries = vectorSpace.getNumGlobalEntries();
        for (int i=0; i < result.length; i++) {
            result[i] /= numGlobalEntries;
        }
        return result;
    }
    
    public void scale(double scalar) {
        Blas myBlas = new NetlibBlas();
        for (int i=0; i < values.length; i++) {
            myBlas.scale(scalar, values[i]);
        }
    }
    
    public double[] dot(MultiVector otherMultiVector) {
        if (!vectorSpace.isCompatible(otherMultiVector.getVectorSpace())) {
            this.println("ERR", "The MultiVectors to be dotted are not compatible.");
        }
        Blas myBlas = new NetlibBlas();
        double[][] y = otherMultiVector.getValues();
        double[] result = new double[vectorSpace.getNumMyGlobalEntries()];
        for(int i=0; i < vectorSpace.getNumMyGlobalEntries(); i++) {
            result[i] = myBlas.dot(this.values[i], y[i]);
        }
        
        result = vectorSpace.getComm().sumAll(result);
        return result;
    }
    
    public double[] norm1() {
        double[] result = new double[values.length];
        Blas myBlas = new NetlibBlas();
        for(int i=0; i < values.length; i++) {
            result[i] = myBlas.asum(values[i]);
        }
        
        result = vectorSpace.getComm().sumAll(result);
        return result;
    }
    
    public double[] norm2() {
        double[] result = new double[values.length];
        int sum;
        for (int i=0; i < values.length; i++) {
            sum = 0;
            for(int j=0; j < values[i].length; j++) {
                sum += values[i][j] * values[i][j];
            }
            result[i] = sum;
        }
        
        result = vectorSpace.getComm().sumAll(result);
        for(int i=0; i < values.length; i++) {
            result[i] = Math.sqrt(result[i]);
        }
        return result;
    }
    
    public double[] normInf() {
        double[] result = new double[values.length];
        
        Blas myBlas = new NetlibBlas();
        int j;
        for (int i=0; i < values.length; i++) {
            j = myBlas.iamax(values[i]);
            result[i] = Math.abs(values[i][j]);
        }
        
        result = vectorSpace.getComm().maxAll(result);
        return result;
    }
    
    public void Reciprocal() {
        for (int i=0; i < values.length; i++) {
            for(int j=0; j < values[i].length; j++) {
                values[i][j] = 1/values[i][j];
            }
        }
    }
    
    public void abs() {
        for (int i=0; i < values.length; i++) {
            for(int j=0; j < values[i].length; j++) {
                values[i][j] = Math.abs(values[i][j]);
            }
        }
    }
    
    public VectorSpace getVectorSpace() {
        return this.vectorSpace;
    }
    
    public double[][] getValues() {
        return this.values;
    }
    
    public void update(double scalarA, MultiVector A, double scalarThis) {
        if (!vectorSpace.isCompatible(A.getVectorSpace())) {
            this.println("ERR", "The MultiVectors are not compatible.");
        }
        
        double[][] valuesA = A.getValues();
        if (scalarThis==0.0) {
            for (int i = 0; i < values.length; i++) {
                for (int j = 0; j < values[i].length; j++) values[i][j] = scalarA * valuesA[i][j];
            }
        }
        else if (scalarThis==1.0) {
            for (int i = 0; i < values.length; i++)
                for (int j = 0; j < values[i].length; j++) {
                    values[i][j] = values[i][j] + scalarA * valuesA[i][j];
                }
        }
        else if (scalarA==1.0) {
            for (int i = 0; i < values.length; i++)
                for (int j = 0; j < values[i].length; j++) {
                    values[i][j] = scalarThis * values[i][j] + valuesA[i][j];
                }
        }
        else {
            for (int i = 0; i < values.length; i++)
                for (int j = 0; j < values[i].length; j++) {
                    values[i][j] = scalarThis * values[i][j] + scalarA *  valuesA[i][j];
                }
        }
    }
    
    public void printOutAllVnodes(String iostream) {
        this.printOut(iostream, true);
    }
    
    public void printOut(String iostream) {
        this.printOut(iostream, false);
    }
    
    public void printOut(String iostream, boolean all) {
        if (all || (vectorSpace.getComm().getVnodeId() == 0)) {
            for(int i=0; i < values.length; i++) {
                for (int j=0; j < values[i].length; j++) {
                    this.println(iostream, "col: " + i + " row: " + j + " value: " + values[i][j]);
                }
                System.out.println("");
            }
        }
    }
    
    public int getNumCols() {
        return this.values.length;
    }
    
    public int getNumRows() {
        if (values.length == 0) {
            return 0;
        }
        return this.values[0].length;
    }
    
    public void copyAndPermute(DistObject distObjectSource, int numSameGids, int[] permuteToLids, int[] permuteFromLids, int combineMode) {
        double[][] srcValues = ((MultiVector) distObjectSource).getValues();
        
        // if srcValues is null, we don't need to do anything
        if (srcValues == null) {
            return;
        }
        
        // if this.values does not exist yet, create it
        if (this.values == null) {
            this.values = new double[srcValues.length][this.vectorSpace.getNumMyGlobalEntries()];
        }
        
        if (combineMode == DistObject.ADD) {
            for(int i=0; i < numSameGids; i++) {
                for(int cols=0; cols < srcValues.length; cols++) {
                    this.values[cols][i] += srcValues[cols][i];
                }
            }
            
            for(int i=0; i < permuteToLids.length; i++) {
                for(int cols=0; cols < srcValues.length; cols++) {
                    this.values[cols][permuteToLids[i]] += srcValues[cols][permuteFromLids[i]];
                }
            }
        }
        else if (combineMode == DistObject.REPLACE) {
            for(int i=0; i < numSameGids; i++) {
                for(int cols=0; cols < srcValues.length; cols++) {
                    this.values[cols][i] = srcValues[cols][i];
                }
            }
            
            for(int i=0; i < permuteToLids.length; i++) {
                for(int cols=0; cols < srcValues.length; cols++) {
                    this.values[cols][permuteToLids[i]] = srcValues[cols][permuteFromLids[i]];
                }
            }
        }
        else {
            this.println("ERR", "The specified combine mode is not supported by MultiVector import/export.");
            System.exit(1);
        }
    }
    
    public Serializable[] packAndPrepare(DistObject distObjectSource, int[] exportGids, int[] exportLids) {
        Serializable[] exportData = new Serializable[exportGids.length];  // the object to be exported by Distributor.distribute
        Serializable[] entryData; // [0] = the Gid, [1] = the values (tmpEntry)
        double[] tmpEntry;  // used to build an array of all the row entries being exported
        double[][] exportValues = ((MultiVector) distObjectSource).getValues();
        for(int i=0; i < exportGids.length; i++) {
            entryData = new Serializable[2];
            exportData[i] = entryData;
            entryData[0] = new Integer(exportGids[i]);  // set the Gid
            tmpEntry = new double[exportValues.length]; // allocate space for all the row entries
            entryData[1] = tmpEntry;  // set the values array
            for(int j=0; j < tmpEntry.length; j++) {
                // copy all the row entries out of exportValues
                tmpEntry[j] = exportValues[j][exportLids[i]];
            }
        }
        
        return exportData;
    }
    
    public int[][] unpackAndCombine(Serializable[] importData, int combineMode) {
        Serializable[] entry;
        Serializable[] entryData;
        int[] intArray;
        int gid;
        int lid;
        double[] importValues;
        
        // if this.values doesn't exist yet, create it
        if (this.values == null) {
            for(int i=0; i < importData.length; i++) {
                // if a vnode didn't send us any data, then importData[vnodeId] == null
                if (importData[i] != null) {
                    entry = (Serializable[]) importData[i];
                    entryData = (Serializable[]) entry[0];
                    importValues = (double[]) entryData[1];
                    int numCols = importValues.length;
                    this.values = new double[numCols][this.vectorSpace.getNumMyGlobalEntries()];
                    break;
                }
            }
            
            // if this.values is still null, then we probably don't have any Gids, or somethin went wrong
            // so assume we don't have any gids
            if (this.values == null) {
                this.values = new double[0][0];
            }
        }
        
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
            reverseExportVnodeIdsGidsLids[0] = new int[vectorSpace.getComm().getNumVnodes()];
            reverseExportVnodeIdsGidsLids[1] = new int[sumEntries];
            reverseExportVnodeIdsGidsLids[2] = new int[sumEntries];
            reverseExportVnodeIdsGidsLids[3] = new int[sumEntries];
        }
        int revCount = 0;
        for(int i=0; i < importData.length; i++) {
            // if a vnode didn't send us any data, then importData[vnodeId] == null
            if (importData[i] != null) {
                entry = (Serializable[]) importData[i];  // get the array of elements sent to us by the vnode i
                // unpack each element
                for(int j=0; j < entry.length; j++) {
                    entryData = (Serializable[]) entry[j];
                    gid = ((Integer) entryData[0]).intValue();
                    lid = vectorSpace.getLocalIndex(gid);
                    
                    //this.println("STD", "MultiVector: Combining Gid: " + gid + " lid: " + lid);
                    // debugging code!!!
                    if (lid == -1) {
                        this.println("ERR", "In MultiVector: Got a bad gid: " + gid);
                        continue;
                    }
                    
                    if (!this.doneForward) {
                        reverseExportVnodeIdsGidsLids[0][i]++;
                        reverseExportVnodeIdsGidsLids[1][revCount] = i;
                        reverseExportVnodeIdsGidsLids[2][revCount] = gid;
                        reverseExportVnodeIdsGidsLids[3][revCount++] = lid;
                    }
                    importValues = (double[]) entryData[1];
                    
                    // combine the existing local values and imported values as defined by the combine mode passed in
                    if (combineMode == DistObject.ADD) {
                        for(int k=0; k < this.values.length; k++) {
                            this.values[k][lid] += importValues[k];
                        }
                    }
                    else if (combineMode == DistObject.REPLACE) {
                        for(int k=0; k < this.values.length; k++) {
                            this.values[k][lid] = importValues[k];
                        }
                    }
                    else {
                        this.println("ERR", "The combine mode you specified is not supported by MultiVector import/export.");
                        System.exit(1);
                    }
                }
            }
        }
        
        this.doneForward = true;
        
        return reverseExportVnodeIdsGidsLids;
    }
    
    public void writeExternal(ObjectOutput out) throws java.io.IOException {
        out.writeObject(this.values);
        out.writeObject(this.vectorSpace);
    }
    
    public void readExternal(ObjectInput in) throws java.io.IOException, ClassNotFoundException {
        this.values = (double[][]) in.readObject();
        this.vectorSpace = (VectorSpace) in.readObject();
    }
    
    public static MultiVector readFromFile(String fileName, Comm comm) {
        try {
            ObjectInputStream ois = new ObjectInputStream(new FileInputStream(fileName));
            MultiVector multiVector = (MultiVector) ois.readObject();
            ois.close();
            multiVector.getVectorSpace().setComm(comm);
            return multiVector;
        } catch (java.io.IOException e) {
            JpetraObject.println("FATALERR", e.toString());
            System.exit(1);
        } catch (java.lang.ClassNotFoundException e) {
            JpetraObject.println("FATALERR", e.toString());
            System.exit(1);
        }
        
        return null; // should never get here
    }
    
    public static void writeToFile(String fileName, MultiVector multiVector) {
        if (multiVector.getVectorSpace().getComm().getVnodeId() != 0) {
            return;  // only vnode 0 will write the object to afile
        }
        if (multiVector.getVectorSpace().getNumMyGlobalEntries() != multiVector.getVectorSpace().getNumGlobalEntries()) {
            JpetraObject.println("FATALERR", "Only a MultiVector with all Gids present on the root vnode (vnode 0) can be written to a serialized object file.");
            System.exit(1);
        }
        
        try {
            ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(fileName));
            oos.writeObject(multiVector);
            oos.close();
        } catch (java.io.IOException e) {
            JpetraObject.println("FATALERR", e.toString());
            System.exit(1);
        }
    }
    
    public Object clone() {
        double[][] cloneValues = new double[this.values.length][];
        for(int i=0; i < this.values.length; i++) {
            cloneValues[i] = new double[this.values[i].length];
            System.arraycopy(this.values[i], 0, cloneValues, 0, this.values[i].length);
        }
        MultiVector cloneMultiVector = new MultiVector(this.vectorSpace, cloneValues);
        cloneMultiVector.doneForward = this.doneForward;
        
        return cloneMultiVector;
    }
    
    public boolean equals(Object obj) {
        // do a direct reference comparison
        if (obj == this) {
            return true;
        }
        MultiVector otherMultiVector = (MultiVector) obj;
        // do a quick length check on the values 2d array
        if (otherMultiVector.getValues().length != this.values.length) {
            this.println("STD", "In MultiVector.equals: otherMultiVector.getValues().length != this.values.length");
            return false;
        }
        // check to see if the vectorSpaces are equal
        if (!otherMultiVector.getVectorSpace().equals(this.vectorSpace)) {
            this.println("STD", "In MultiVector.equals: !otherMultiVector.getVectorSpace().equals(this.vectorSpace)");
            return false;
        }
        // check to see if the values 2d arrays are the same
        
        /*double[][] otherValues = otherMultiVector.getValues();
        for(int i=0; i < this.values.length; i++) {
            if (this.values[i].length != otherValues[i].length) {
                return false;
            }
            for(int j=0; j < this.values[i].length; j++) {
                if (this.values[i][j] != otherValues[i][j]) {
                    return false;
                }
            }
        }*/
        
        double[][] otherValues = otherMultiVector.getValues();
        for(int i=0; i < this.values.length; i++) {
            if (!java.util.Arrays.equals(this.values[i], otherValues[i])) {
                this.println("STD", "In MultiVector.equals: !java.util.Arrays.equals(this.values[i], otherValues[i])");
                return false;
            }
        }
        
        return true;
    }
}
