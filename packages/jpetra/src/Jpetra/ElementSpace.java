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

import java.util.TreeMap;

/**
 *
 * @author  Jason Cross
 */
public class ElementSpace {
    private int[] myGlobalElements;
    private Comm comm;
    private TreeMap gidsToLids;
    private Directory directory;
    
    // values that describe the distribution of Global and Local Elements
    private int numGlobalElements;
    private int numMyGlobalElements;
    private int minMyGlobalElementId;
    private int maxMyGlobalElementId;
    private int minLocalElementId;
    private int maxLocalElementId;
    private int minGlobalElementId;
    private int maxGlobalElementId;
    
    // values that are needed to compute vnodeIds from Gids
    private int numRemainderIndices;
    private int numIndicesPerVnode;
    
    private boolean distributedGlobally;
    private boolean distributedLinearly;
    /**
     * Constructs an <code>ElementSpace</code> automatically based on the contigous
     * even distribution of global elements to each vnode.
     */
    public ElementSpace(int numGlobalElements, int indexBase, Comm comm) {
        // for both Serial and Parallel
        this.comm = comm;
        this.numGlobalElements = numGlobalElements;
        
        // for serial only
        if (comm.isSerial()) {
            // do nothing
            this.minMyGlobalElementId = indexBase;
            this.numMyGlobalElements = this.numGlobalElements;
            this.distributedGlobally = false;
            this.distributedLinearly = true;
            return;
        }
        
        // for parallel only
        // !! right now this is for linear uniform distributions!
        this.distributedGlobally = true;
        this.distributedLinearly = true;
        numIndicesPerVnode = numGlobalElements / comm.getNumVnodes();
        numRemainderIndices = numGlobalElements % comm.getNumVnodes();
        this.numMyGlobalElements = numIndicesPerVnode;
        int numRemainderIndicesToAdd = 0;
        // if the ElementSpace does not map evenly onto all vectors
        // then give vnodes < remainder an additional element
        if (comm.getVnodeId() < this.numRemainderIndices) {
            this.numMyGlobalElements++;
            // accounts for the indices owned by vnodes < myVnode
            numRemainderIndicesToAdd = comm.getVnodeId();
        }
        else {
            numRemainderIndicesToAdd = numRemainderIndices;
        }
        this.minMyGlobalElementId = (comm.getVnodeId() * this.numIndicesPerVnode) + this.numRemainderIndices + indexBase;
    }
    
    /**
     * Constructs an arbitrarly defined <code>ElementSpace</code> using <code>myGlobalElements</code>
     * passed in by the user.
     */
    public ElementSpace(int[] myGlobalElements, Comm comm) {
        // for both serial and parallel
        this.distributedLinearly = false;
        this.myGlobalElements = myGlobalElements;
        this.numMyGlobalElements = myGlobalElements.length;
        this.minLocalElementId = 0;
        this.maxLocalElementId = myGlobalElements.length - 1;
        
        // build Gids To Lids Mapping
        this.gidsToLids = new TreeMap();
        for(int i=0; i < myGlobalElements.length; i++) {
            gidsToLids.put(new Integer(myGlobalElements[i]), new Integer(i));
        }
        
        // find my min and max global elements
        int max = myGlobalElements[0];
        int min = myGlobalElements[0];
        for(int i=1; i < myGlobalElements.length; i++) {
            if (max < myGlobalElements[i]) {
                max = myGlobalElements[i];
            }
            if (min > myGlobalElements[i]) {
                min = myGlobalElements[i];
            }
        }
        this.minMyGlobalElementId = min;
        this.maxMyGlobalElementId = max;
        this.comm = comm;
        
        // serial only
        if (comm.isSerial()) {
            this.numGlobalElements = myGlobalElements.length;
            this.minGlobalElementId = this.minMyGlobalElementId;
            this.maxGlobalElementId = this.maxMyGlobalElementId;
            this.distributedGlobally = false;
            return;
        }
        
        // parallel only
        this.distributedGlobally = true;
        // need to get the number of global elements by doing a gather of the number
        // of elements on each vnode
        int[] tmp = comm.sumAll(new int[]{this.numMyGlobalElements});
        this.numGlobalElements = tmp[0];
        // find global min and max for global elements
        tmp = new int[]{this.minMyGlobalElementId};
        tmp = comm.minAll(tmp);
        this.minGlobalElementId = tmp[0];
        tmp = new int[]{this.maxMyGlobalElementId};
        tmp = comm.maxAll(tmp);
        this.maxGlobalElementId = tmp[0];
    }
    
    /**
     * Constructs an <code>ElementSpace</code> automatically based on the user defined
     * linear distribution of global elements given to each vnode.
     *
     * @param numGlobalElements -1 for globaly distributed or numMyGlobalElements for local/serial
     */
    public ElementSpace(int numGlobalElements, int numMyGlobalElements, int indexBase, Comm comm) {
        this.numMyGlobalElements = numMyGlobalElements;
        this.comm = comm;
        this.distributedLinearly = true;
        
        if (comm.isSerial() || (numGlobalElements == numMyGlobalElements)) {
            this.distributedGlobally = false;
            this.numGlobalElements = numGlobalElements;
            this.numMyGlobalElements = numMyGlobalElements;
            this.minMyGlobalElementId = indexBase;
            this.minGlobalElementId = this.minMyGlobalElementId;
            this.maxMyGlobalElementId = numGlobalElements + indexBase;
            this.maxGlobalElementId = this.maxMyGlobalElementId;
            return;
        }
        
        // find the number of global elements by finding the sum of all local elements
        int[] sum = comm.sumAll(new int[]{numMyGlobalElements});
        this.numGlobalElements = sum[0];
        int[] scanSum = comm.scanSums(new int[]{numMyGlobalElements});
        int myMaxGlobalElementId = scanSum[0] + indexBase;
        this.minMyGlobalElementId = myMaxGlobalElementId - numMyGlobalElements + indexBase;
        
    }
    
    
    public int getLocalElementId(int globalElementId) {
        // for both serial and parallel
        if ((globalElementId < this.minGlobalElementId) || (globalElementId > this.maxGlobalElementId)) {
            return -1;
        }
        
        // !! need to add support for arbitrarly defined element spaces
        
        // for even contigous distribution compute the global ID
        // serial only
        if (comm.isSerial()) {
            return globalElementId - minMyGlobalElementId;
        }
        
        // linear continous parallel only
        if (this.distributedLinearly) {
            return globalElementId - minMyGlobalElementId;
        }
        
        // abritrarily defined parrallel
        Object result = gidsToLids.get(new Integer(globalElementId));
        if (result == null) {
            return -1;
        }
        return ((Integer) result).intValue();
    }
    
    public int getGlobalElementId(int localElementId) {
        // for both serial and parallel
        if ((localElementId < this.minLocalElementId) || (localElementId > this.maxLocalElementId)) {
            return -1;
        }
        
        // !! need to add support for arbitrarly defined element spaces
        
        // for even continous distribution compute the global ID
        // serial only
        if (comm.isSerial()) {
            return localElementId + minMyGlobalElementId;
        }
        
        // liner continous parallel only
        if (this.distributedLinearly) {
            return minMyGlobalElementId + localElementId;
        }
        
        // abritrarily defined parrallel
        return this.myGlobalElements[localElementId];
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
    
    public int[] getMyGlobalElementIds() {
        if (this.myGlobalElements != null) {
            return this.myGlobalElements;
        }
        
        int[] generatedGlobalElementIds = null;
        if (distributedLinearly) {
            generatedGlobalElementIds = new int[this.numMyGlobalElements];
            for (int i=0; i < this.numMyGlobalElements; i++) {
                generatedGlobalElementIds[i] = i + this.minMyGlobalElementId;
            }
        }
        
        return generatedGlobalElementIds;
    }
    
    /*public int[] getRemoteVnodeIdList(int[] remoteGlobalElementIds, int[] remoteVnodeIds, int[] remoteLocalElementIds) {
        if (this.directory == null) {
            this.directory = new BasicDirectory();
        }
        return null; // temporary so class will compile
    }*/
    
    public boolean isDistributedGlobally() {
        return this.distributedGlobally;
    }
    
    public boolean isDistributedLinearly() {
        return this.distributedLinearly;
    }
    
    public Comm getComm() {
        return this.comm;
    }
    
    public int getNumRemainderIndices() {
        return this.numRemainderIndices;
    }
    
    public int getNumIndicesPerVnode() {
        return this.numIndicesPerVnode;
    }
}
