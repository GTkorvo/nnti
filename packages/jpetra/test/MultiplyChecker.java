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
        // setup an output stream for main output
        this.outputStreams.put("MAIN", new Output("-->", true, System.out, false, System.out));
        
        //Comm comm = new CcjComm("test/ccjhosts.txt");
        Comm comm = new SerialComm();
        
        // read in in matrix A and distribute it
        CisMatrix inMatrix = null;
        try {
            this.println("MAIN", "reading in inMatrix from " + matrixAfileName + "...");
            inMatrix = CisMatrixReader.read(matrixAfileName, rowOriented, comm);
        } catch (java.io.IOException e) {
            this.println("FATALERR", e.toString());
            System.exit(1);
        }
        
        this.println("MAIN", "inMatrix.getPrimaryVectorSpace().getNumGlobalEntries(): " + inMatrix.getPrimaryVectorSpace().getNumGlobalEntries());
        VectorSpace globalVectorSpace = new VectorSpace(new ElementSpace(inMatrix.getPrimaryVectorSpace().getNumGlobalEntries(), 0, comm));
        CisMatrix globalMatrix = new CisMatrix(globalVectorSpace, CisMatrix.ROW_ORIENTED);
        Export exporter = new Export(inMatrix.getPrimaryVectorSpace(), globalVectorSpace);
        globalMatrix.exportValues(inMatrix, exporter, DistObject.REPLACE);
        globalMatrix.fillComplete();
        
        // readin in multivector x and distribute it
        MultiVector xInMultiVector = null;
        try {
            this.println("MAIN", "reading in multivector x from " + multivectorXfileName + "...");
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
        this.println("MAIN", "printing result MultiVector y");
        y.printOutAllVnodes("MAIN");
        
        this.println("MAIN", "Multiplication y=x*A is complete.  Now import all values from y into outMatrix on vnode 0...");
        
        MultiVector outMatrix;
        if (comm.getVnodeId() == 0) {
            double[][] outMatrixValues = new double[globalX.getNumCols()][globalMatrix.getRowVectorSpace().getNumGlobalEntries()];
            outMatrix = new MultiVector(new VectorSpace(new ElementSpace(y.getVectorSpace().getNumGlobalEntries(), y.getVectorSpace().getNumGlobalEntries(), 0, comm)), outMatrixValues);
        }
        else {
            outMatrix = new MultiVector(new VectorSpace(new ElementSpace(y.getVectorSpace().getNumGlobalEntries(), 0, 0, comm)), new double[0][0]);
        }
        this.println("MAIN", "outMatrixNumRows: " + outMatrix.getNumRows() + " outMatrixNumMyGlobalElementIds: " + outMatrix.getVectorSpace().getNumMyGlobalEntries());
        Import importer = new Import(y.getVectorSpace(), outMatrix.getVectorSpace());
        outMatrix.importValues(y, importer, DistObject.REPLACE);
        this.println("MAIN", "Import of y to outMatrix on vnode 0 is complete.");
        
        if (comm.getVnodeId() == 0) {
            this.println("MAIN", "Checking result with multivector B...");
        }
        
        MultiVector bMutliVector = null;
        try {
            this.println("MAIN", "reading in multivector b from " + multivectorBfileName + "...");
            bMutliVector = MultiVectorReader.read(multivectorBfileName, comm);
        } catch (java.io.IOException e) {
            this.println("FATALERR", e.toString());
            System.exit(1);
        }
        
        /*
        try {
            MultiVectorWriter.write("javaout", outMatrix, comm);
        }
        catch (java.io.IOException e) {
            this.println("ERR", e.toString());
        }
        */
        
        if (comm.getVnodeId() == 0) {
            if (bMutliVector.equals(outMatrix)) {
                this.println("MAIN", "Test Passed!");
            }
            else {
                this.println("MAIN", "Test Failed!");
                double[][] computedBvalues = outMatrix.getValues();
                double[][] bValues = bMutliVector.getValues();
                if (computedBvalues.length != bValues.length) {
                    this.println("FATALERR", "The computed b and the read in b have different number of rows.");
                    System.exit(1);
                }
                else if (computedBvalues[0].length != bValues[0].length) {
                    this.println("FATALERR", "The computed b and the read in b have different number of columns.");
                    System.exit(1);
                }
                
                String diff;
                String noDiff = "";
                String differ = " VALUES DIFFER!";
                for(int i=0; i < computedBvalues.length; i++) {
                    for(int j=0; j < computedBvalues[i].length; j++) {
                        if (computedBvalues[i][j] == bValues[i][j]) {
                            diff = noDiff;
                        }
                        else {
                            diff = differ;
                        }
                        this.println("MAIN", "i: " + i + " j: " + j + " b: " + bValues[i][j] + " computB: " + computedBvalues[i][j] + diff);
                    }
                }
               /* this.println("MAIN", "printing result matrix...");
                outMatrix.printOut("MAIN");
                this.println("MAIN", "printing b...");
                bMutliVector.printOut("MAIN"); */
            }
        }
        
        
        this.println("MAIN", "Done!");
        System.exit(0);
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        if (args.length != 5) {
            System.out.println("Usage:");
            System.out.println("test/MultiplyChecker matrixA.mtx multivectorX.mtx multivectorB.mtx [row|column] [no_transpose|transpose]");
            System.exit(1);
        }
        
        boolean rowOriented = CisMatrix.COL_ORIENTED;
        if (args[3].equals("row")) {
            rowOriented = CisMatrix.ROW_ORIENTED;
        }
        else if (!args[3].equals("column")) {
            System.out.println("Error: You must specify whether matrix A is row or column oriented.");
            System.exit(1);
        }
        
        boolean useTranspose = CisMatrix.USE_A;
        if (args[4].equals("transpose")) {
            useTranspose = CisMatrix.USE_TRANSPOSE_A;
        }
        else if (!args[4].equals("no_transpose")) {
            System.out.println("Error: You must specify whether to use the tranpose of matrix A or not.");
            System.exit(1);
        }
        
        new MultiplyChecker(args[0], args[1], args[2], rowOriented, useTranspose);
    }
    
}
