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

package test;

import Jpetra.*;
import Jpetra.MatrixMarketIO.*;

/**
 *
 * @author  Jason Cross
 */
public class MultiplyChecker extends JpetraObject {
    public MultiplyChecker(String matrixAfileName, String multivectorXfileName, String multivectorBfileName, boolean rowOriented, boolean useTranspose) {
        initializeOutput();
        Comm comm = new CcjComm("test/ccjhosts.txt");
        
        // read in in matrix A and distribute it
        CisMatrix inMatrix = null;
        try {
            this.println("STD", "reading in inMatrix from " + matrixAfileName + "...");
            inMatrix = CisMatrixReader.read(matrixAfileName, rowOriented, comm);
        } catch (java.io.IOException e) {
            this.println("FATALERR", e.toString());
            System.exit(1);
        }
        
        VectorSpace globalVectorSpace = new VectorSpace(new ElementSpace(inMatrix.getPrimaryVectorSpace().getNumGlobalEntries(), 0, comm));
        CisMatrix globalMatrix = new CisMatrix(globalVectorSpace, CisMatrix.ROW_ORIENTED);
        Export exporter = new Export(inMatrix.getPrimaryVectorSpace(), globalVectorSpace);
        globalMatrix.exportValues(inMatrix, exporter, DistObject.REPLACE);
        globalMatrix.fillComplete();
        
        // readin in multivector x and distribute it
        MultiVector xInMultiVector = null;
        try {
            this.println("STD", "reading in multivector x from " + multivectorXfileName + "...");
            xInMultiVector = MultiVectorReader.read(multivectorXfileName, comm);
        } catch (java.io.IOException e) {
            this.println("FATALERR", e.toString());
            System.exit(1);
        }
        
        VectorSpace globalXvs = new VectorSpace(new ElementSpace(xInMultiVector.getVectorSpace().getNumGlobalEntries(), 0, comm));
        MultiVector globalX = new MultiVector(globalXvs);
        Export exporter2 = new Export(xInMultiVector.getVectorSpace(), globalXvs);
        globalX.exportValues(xInMultiVector, exporter2, DistObject.REPLACE);
        
        double[][] yValues = new double[globalX.getNumCols()][globalVectorSpace.getNumMyGlobalEntries()];
        MultiVector y = new MultiVector(globalVectorSpace, yValues);
        
        globalMatrix.multiply(useTranspose, globalX, y);
        this.println("STD", "-->printing result MultiVector y");
        y.printOutAllVnodes("STD");
        
        this.println("STD", "-->Multiplication y=x*A is complete.  Now import all values from y into outMatrix on vnode 0...");
        
        MultiVector outMatrix;
        if (comm.getVnodeId() == 0) {
            double[][] outMatrixValues = new double[globalX.getNumCols()][globalMatrix.getRowVectorSpace().getNumGlobalEntries()];
            outMatrix = new MultiVector(new VectorSpace(new ElementSpace(y.getVectorSpace().getNumGlobalEntries(), y.getVectorSpace().getNumGlobalEntries(), 0, comm)), outMatrixValues);
        }
        else {
            outMatrix = new MultiVector(new VectorSpace(new ElementSpace(y.getVectorSpace().getNumGlobalEntries(), 0, 0, comm)), new double[0][0]);
        }
        this.println("STD", "-->outMatrixNumRows: " + outMatrix.getNumRows() + " outMatrixNumMyGlobalElementIds: " + outMatrix.getVectorSpace().getNumMyGlobalEntries());
        Import importer = new Import(y.getVectorSpace(), outMatrix.getVectorSpace());
        outMatrix.importValues(y, importer, DistObject.REPLACE);
        this.println("STD", "-->Import of y to outMatrix on vnode 0 is complete.");
        
        if (comm.getVnodeId() == 0) {
            this.println("STD", "-->Checking result with multivector B.");
            MultiVector bMutliVector = null;
            try {
                this.println("STD", "reading in multivector b from " + multivectorBfileName + "...");
                bMutliVector = MultiVectorReader.read(multivectorBfileName, comm);
            } catch (java.io.IOException e) {
                this.println("FATALERR", e.toString());
                System.exit(1);
            }
            
            if (bMutliVector.equals(outMatrix)) {
                this.println("STD", "-->Test Passed!");
            }
            else {
                this.println("STD", "--Test Failed!");
                this.println("STD", "printing result matrix...");
                outMatrix.printOut("STD");
                this.println("STD", "printing b...");
                bMutliVector.printOut("STD");
            }
        }
        
        
        this.println("STD", "-->Done!");
        System.exit(0);
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        if (args.length != 5) {
            System.out.println("Usage:");
            System.out.println("test/MultiplyChecker matrixA.mtx multivectorX.mtx multivectorB.mtx [row|column] [no_transpose|transpose]");
        }
        boolean rowOriented = CisMatrix.COL_ORIENTED;
        if (args[3].equals("row")) {
            rowOriented = CisMatrix.ROW_ORIENTED;
        }
        else if (!args[3].equals("column")) {
            System.out.println("Error: You must specify whether matrix A is row or column oriented.");
        }
        
        boolean useTranspose = CisMatrix.USE_A;
        if (args[4].equals("transpose")) {
            useTranspose = CisMatrix.USE_TRANSPOSE_A;
        }
        else if (!args[4].equals("no_transpose")) {
            System.out.println("Error: You must specify whether to use the tranpose of matrix A or not.");
        }
        
        new MultiplyChecker(args[0], args[1], args[2], rowOriented, useTranspose);
    }
    
}
