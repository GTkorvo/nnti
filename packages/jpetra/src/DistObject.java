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

public class DistObject extends JpetraObject {
   // int err;
    
    int numRemoteIDs;
    
    int numExportIDs;
    int[] exportGIDs;
    int[] exportPIDs;
    
    public DistObject (int numRemoteIDs) {
        //this.err  = err;
        this.numRemoteIDs = numRemoteIDs;   
    }

    public DistObject (int numExportIDs, int[] exportGIDs, int[] exportPIDs) {
        //this.err  = err;
        this.numExportIDs = numExportIDs;
        this.exportGIDs = exportGIDs;
        this.exportPIDs = exportPIDs;
    }
    
   // public int getErr () {
   //     return this.err;
   // }
    
    public int getNumExportIDs() {
        return numExportIDs;
    }
    
        public int[] getExportGIDs() {
        return exportGIDs;
    }
    
        public int[] getExportPIDs() {
        return exportPIDs;
    }
    
        public int getNumRemoteIDs() {
        return numRemoteIDs;
    }
    
}
