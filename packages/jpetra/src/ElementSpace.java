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

/**
 *
 * @author  Jason Cross
 */
public class ElementSpace {
    private int[] myGlobalElements;
    private Comm comm;
    
    // values that describe the distribution of Global and Local Elements
    private int numGlobalElements;
    private int numMyGlobalElements;
    private int numIndicesPerVnode;
    private int numRemainderIndices;
    private int myGlobalStartIndex;
    private int minLocalElementId;
    private int maxLocalElementId;
    private int minGlobalElementId;
    private int maxGlobalElementId;
    
    /**
     * Constructs an <code>ElementSpace</code> automatically based on the contigous
     * even distribution of global elements to each vnode.
     */
    public ElementSpace(int numGlobalElements, Comm comm) {
        // for both Serial and Parallel
        this.comm = comm;
        this.numGlobalElements = numGlobalElements;
        
        // for serial only
        if (comm.getIsSerial()) {
            // do nothing
            this.numMyGlobalElements = this.numGlobalElements;
            return;
        }
        
        // for parallel only
        numIndicesPerVnode = numGlobalElements / comm.getNumVnodes();
        numRemainderIndices = numGlobalElements % comm.getNumVnodes();
        this.numMyGlobalElements = numIndicesPerVnode;
        int numRemainderIndicesToAdd = 0;
        // if the ElementSpace does not map evenly onto all vectors
        // then give vnodes < remainder an additional element
        if (comm.getVnodeID() < this.numRemainderIndices) {
            this.numMyGlobalElements++;
            // accounts for the indices owned by vnodes < myVnode
            numRemainderIndicesToAdd = comm.getVnodeID();
        }
        else {
            numRemainderIndicesToAdd = numRemainderIndices;
        }
        this.myGlobalStartIndex = (comm.getVnodeID() * this.numIndicesPerVnode) + this.numRemainderIndices;
    }
    
    /**
     * Constructs an arbitrarly defined <code>ElementSpace</code> using <code>myGlobalElements</code>
     * passed in by the user.
     */
    public ElementSpace(int[] myGlobalElements, Comm comm) {
        // for both serial and parallel
        this.myGlobalElements = myGlobalElements;
        this.numMyGlobalElements = myGlobalElements.length;
        this.comm = comm;
        
        // serial only
        if (comm.getIsSerial()) {
            this.numGlobalElements = myGlobalElements.length;
            return;
        }
        
        // parallel only
        // need to get the number of global elements by doing a gather of the number
        // of elements on each vnode
        // this.numGlobalElements = ?
    }
    
    public int getLocalElementId(int globalElementId) {
        // for both serial and parallel
        if ((globalElementId < minGlobalElementId) || (globalElementId > maxGlobalElementId)) {
            return -1;
        }
        
        // !! need to add support for arbitrarly defined element spaces
        
        // for even contigous distribution compute the global ID
        // serial only
        if (comm.getIsSerial()) {
            return globalElementId;
        }
        
        // parallel only
        return globalElementId - myGlobalStartIndex;
    }
    
    public int getGlobalElementId(int localElementId) {
        // for both serial and parallel
        if ((localElementId < minLocalElementId) || (localElementId > maxLocalElementId)) {
            return -1;
        }
        
        // !! need to add support for arbitrarly defined element spaces
        
        // for even contigous distribution compute the global ID
        // serial only
        if (comm.getIsSerial()) {
            return localElementId;
        }
        
        // parallel only
        return myGlobalStartIndex + localElementId;
    }
    
    public boolean isMyLocalElementId(int globalElementId) {
        if (this.getLocalElementId(globalElementId) != -1) {
            return true;
        }
        return false;
    }
    
    public boolean isMyGlobalElementId(int localElementId) {
        if (this.getLocalElementId(localElementId) != -1) {
            return true;
        }
        return false;
    }
    
    public int getNumGlobalElements() {
        return this.numGlobalElements;
    }
    
    public int getNumMyGlobalElements() {
        return this.numMyGlobalElements;
    }
    
    public int getMinLocalElementId() {
        return this.minLocalElementId;
    }
    
    public int getMaxLocalElementId() {
        return this.maxLocalElementId;
    }
    
    public int getMinGlobalElementId() {
        return this.minGlobalElementId;
    }
    
    public int getMaxGlobalElementId() {
        return this.maxGlobalElementId;
    }
    
    public Comm getComm() {
        return this.comm;
    }
}
