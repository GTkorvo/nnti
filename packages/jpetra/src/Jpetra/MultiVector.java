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

/**
 *
 * @author  Jason Cross
 */
public class MultiVector extends DistObject {
    VectorSpace vectorSpace;
    double[][] values;  // in column major form
    
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
    
    public void printOut(String iostream) {
        for(int i=0; i < values.length; i++) {
            for (int j=0; j < values[i].length; j++) {
                this.println(iostream, "col: " + i + " row: " + j + " value: " + values[i][j]);
            }
            System.out.println("");
        }
    }
    
    public int getNumCols() {
        return this.values.length;
    }
    
    public int getNumRows() {
        return this.values[0].length;
    }
    
    public void copyAndPermute(DistObject distObjectSource, int numSameGids, int[] permuteToLids, int[] permuteFromLids) {
        double[][] srcValues = ((MultiVector) distObjectSource).getValues();
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
    
    public void unpackAndCombine(Serializable[] importData, int combineMode) {
        Serializable[] entry;
        Serializable[] entryData;
        int[] intArray;
        int gid;
        int lid;
        double[] importValues;
        for(int i=0; i < importData.length; i++) {
            // if a vnode didn't send us any data, then importData[vnodeId] == null
            if (importData[i] != null) {
                entry = (Serializable[]) importData[i];  // get the array of elements send to us by the vnode i
                // unpack each element
                for(int j=0; j < entry.length; j++) {
                    entryData = (Serializable[]) entry[j];
                    gid = ((Integer) entryData[0]).intValue();
                    lid = vectorSpace.getLocalIndex(gid);
                    importValues = (double[]) entryData[1];
                    
                    // combine the existing local values and imported values as defined by the combine mode passed in
                    if (combineMode == DistObject.ADD) {
                        for(int k=0; k < this.values.length; k++) {
                            this.values[k][lid] += importValues[k];
                        }
                    }
                    else {
                        this.println("ERR", "The combine mode you specified is not support yet!");
                        System.exit(0);
                    }
                }
            }
        }
    }
}
