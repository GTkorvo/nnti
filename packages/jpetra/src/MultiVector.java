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

import com.braju.format.*;
/**
 *
 * @author  Jason Cross
 */
public class MultiVector extends DistObject {
    VectorSpace vectorSpace;
    double[][] values;
    
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
                //this.println(iostream, i + " " + j + " " + values[i][j]);
                Format.printf("%' 'f ", new Parameters(values[i][j]));
            }
            System.out.println("");
        }
    }
}
