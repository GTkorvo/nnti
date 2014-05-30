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

import java.io.Serializable;

/**
 *
 * @author  Jason Cross
 */
public interface Distributor {
    /**
     *
     * @param remoteGlobalElementIds Array of GlobalElementIds that this vnode wants.
     * @param remoteGlobalVnodeIds Array of vnode IDs that correspond to the vnodes that own the global elements specified by <code>remoteGlobalElementIds</code>.
     * @return Number of global elements this vnode will receive.
     */
    //public int[] createFromRecieves(int[] remoteGlobalElementIds, int[] remoteGlobalVnodeIds, int[] exportElementIds, int[] exportVnodeIds);
    public int[] createFromReceives(int[] remoteGids, int[] remoteVnodeIds, Comm comm);
    
    /**
     *
     * @param exportVnodeIds The vnodes to export my global elements to.
     */
    public void createFromSends(int[] exportVnodeIds, Comm comm);
    
    public Serializable[] distribute(Serializable[] exportObjects, boolean doReverse);
    public int[][] distribute(int[] toSendData);
    
    public int[] getSenders();
    public int[] getExportVnodeIds();
    
    
    // for reverse op

    public void setReverseExportVnodeIdsGidsLids(int[][] reverseExportVnodeIdsGidsLids);
    
    public int[] getReverseExportVnodeIds();
    
    public int[] getReverseExportGids();
    
    public int[] getReverseExportLids();
    
    public boolean doneForwardOp();
    
    public void setDoneForwardOp(boolean doneForwardOp);
}
