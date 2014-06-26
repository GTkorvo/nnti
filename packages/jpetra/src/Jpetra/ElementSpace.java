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

import java.util.TreeMap;
import java.io.Externalizable;
import java.io.ObjectOutput;
import java.io.ObjectInput;
import java.io.Serializable;

/**
 *
 * @author  Jason Cross
 */
public class ElementSpace extends JpetraObject implements Externalizable {
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
    private boolean distributedUniformly;
    
    public ElementSpace() {
        //empty
    }
    
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
            this.minMyGlobalElementId = indexBase;
            this.numMyGlobalElements = this.numGlobalElements;
            this.maxMyGlobalElementId = this.minMyGlobalElementId + this.numMyGlobalElements-1;
            this.minGlobalElementId = this.minMyGlobalElementId;
            this.maxGlobalElementId = this.maxMyGlobalElementId;
            this.distributedGlobally = false;
            this.distributedLinearly = true;
            this.distributedUniformly = true;
            return;
        }
        
        // for parallel only
        // !! right now this is for linear uniform distributions!
        this.distributedGlobally = true;
        this.distributedLinearly = true;
        this.distributedUniformly = true;
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
        this.minMyGlobalElementId = (comm.getVnodeId() * this.numIndicesPerVnode) + numRemainderIndicesToAdd + indexBase;
        this.maxMyGlobalElementId = this.minMyGlobalElementId + this.numMyGlobalElements-1;
        this.minGlobalElementId = indexBase;
        this.maxGlobalElementId = this.minGlobalElementId + this.numGlobalElements-1;
    }
    
    /**
     * Constructs an arbitrarily defined <code>ElementSpace</code> using <code>myGlobalElements</code>
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
        
        
        Object firstKey = null;
        try {
            firstKey = this.gidsToLids.firstKey();
        }
        catch (java.util.NoSuchElementException e) {
            // empty
        }
        
        if (firstKey != null) {
            this.minMyGlobalElementId = ((Integer) firstKey).intValue();
            this.maxMyGlobalElementId = ((Integer) this.gidsToLids.lastKey()).intValue();
        }
        
        /*
        // should be able to use TreeMap form above to do this!!!!!!!!!!
        // so fix it!!!!
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
         */
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
        
        /* this probably isn't the way to do it, use global min and max to compute total number of elements
        // need to get the number of global elements by doing a gather of the number
        // of elements on each vnode
        int[] tmp = comm.sumAll(new int[]{this.numMyGlobalElements});
        this.numGlobalElements = tmp[0];
         */
        
        
        // find global min and max for global elements
        //int[] tmp = new int[]{this.minMyGlobalElementId};
        //tmp = comm.minAll(tmp);
        //this.minGlobalElementId = tmp[0];
        this.minGlobalElementId = 0; // assume index base is 0
        int[] tmp = new int[]{this.maxMyGlobalElementId};
        tmp = comm.maxAll(tmp);
        this.maxGlobalElementId = tmp[0];
        // find the number of global elements
        //this.println("STD", "ElementSpace: maxGid: " + this.maxGlobalElementId + " minGid: " + this.minGlobalElementId);
        this.numGlobalElements = this.maxGlobalElementId - this.minGlobalElementId + 1;
        
        // find global min and max for global elements
        /*Integer minInteger = null;
        Integer maxInteger = null;
        if (this.numMyGlobalElements != 0) {
            minInteger = new Integer(this.minMyGlobalElementId);
            maxInteger = new Integer(this.maxMyGlobalElementId);
        }*/
        
        //this.println("STD", "ElementSpace(int[], Comm): minMyGlobalElementId: " + this.minMyGlobalElementId + " maxMyGlobalElementId: " + this.maxMyGlobalElementId);
        /*Serializable tmp = comm.minAll((Serializable) minInteger);
        this.minGlobalElementId = ((Integer) tmp).intValue();
        tmp = comm.maxAll((Serializable) maxInteger);
        this.maxGlobalElementId = ((Integer) tmp).intValue();
         */
        // find the number of global elements
        //this.println("STD", "ElementSpace(int[], Comm): maxGid: " + this.maxGlobalElementId + " minGid: " + this.minGlobalElementId);
        this.numGlobalElements = this.maxGlobalElementId - this.minGlobalElementId + 1;
    }
    
    /**
     * Constructs an <code>ElementSpace</code> automatically based on the user defined
     * linear distribution of global elements given to each vnode.
     *
     * @param numGlobalElements -1 to have the Comm determine, actually this parameter does nothing right now...
     *//*
     *                          if numGlobalElements == numMyGlobalElements then 
     *                          the ElementSpace will be serial
     */
    public ElementSpace(int numGlobalElements, int numMyGlobalElements, int indexBase, Comm comm) {
        this.numMyGlobalElements = numMyGlobalElements;
        this.minGlobalElementId = indexBase;
        this.comm = comm;
        this.distributedLinearly = true;
        this.distributedUniformly = false;
         
        //if (comm.isSerial() || (numGlobalElements == numMyGlobalElements)) {
        if (comm.isSerial()) {
            this.distributedGlobally = false;
            this.numGlobalElements = this.numMyGlobalElements;
            this.minMyGlobalElementId = this.minGlobalElementId;
            //this.minGlobalElementId = this.minMyGlobalElementId;
            this.maxMyGlobalElementId = this.numMyGlobalElements + indexBase - 1;
            this.maxGlobalElementId = this.maxMyGlobalElementId;
            return;
        }
        
        // a serial ElementSpace that acts like its global, can be used for multiplication, etc
        if (numGlobalElements == numMyGlobalElements) {
            this.distributedGlobally = true;
            this.numGlobalElements = this.numMyGlobalElements;
            this.minMyGlobalElementId = this.minGlobalElementId;
            //this.minGlobalElementId = this.minMyGlobalElementId;
            this.maxMyGlobalElementId = this.numMyGlobalElements + indexBase - 1;
            this.maxGlobalElementId = this.maxMyGlobalElementId;
            return;
        }
        
        this.distributedGlobally = true;
        // find the number of global elements by finding the sum of all local elements
        int[] sum = comm.sumAll(new int[]{this.numMyGlobalElements});
        this.numGlobalElements = sum[0];
        this.maxGlobalElementId = this.numGlobalElements - 1;
        int[] scanSum = comm.scanSums(new int[]{this.numMyGlobalElements});
        //this.println("STD", "ElementSpace: there are " + scanSum[0] + " gids before mine. I have " + this.numMyGlobalElements + " gids.");
        this.maxMyGlobalElementId = scanSum[0] + indexBase - 1;  // -1 to make indexing start at 0
        this.minMyGlobalElementId = this.maxMyGlobalElementId - this.numMyGlobalElements + indexBase;
        if (comm.getVnodeId() == 0) {
            this.minMyGlobalElementId = indexBase;
        }
        //this.println("STD", "In ElementSpace: myMaxGlobalElementId: " + this.maxMyGlobalElementId + " minMyGlobalElementId: " + this.minMyGlobalElementId + " numMyGids: " + this.numMyGlobalElements + " numGids: " + this.numGlobalElements);
    }
    
    
    public int getLocalElementId(int globalElementId) {
        // for both serial and parallel
        if(this.numMyGlobalElements == 0) {
            return -1;
        }
        
        if ((globalElementId < this.minMyGlobalElementId) || (globalElementId > this.maxMyGlobalElementId)) {
            //this.println("STD", "GlobalId " + globalElementId + " failed min/max (" + this.minMyGlobalElementId + "/" + this.maxMyGlobalElementId + ") test. getNumMyGlobalElements(): " + this.getNumMyGlobalElements());
            return -1;
        }
        
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
            //this.println("STD", "GlobalId " + globalElementId + " failed to be found in gidsToLids set.");
            return -1;
        }
        return ((Integer) result).intValue();
    }
    
    public int getGlobalElementId(int localElementId) {
        // for both serial and parallel
        if ((localElementId < this.minLocalElementId) || (localElementId > this.maxLocalElementId)) {
            return -1;
        }
        
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
    
    public int getMyMinGlobalElementId() {
        return this.minMyGlobalElementId;
    }
    
    public int getMyMaxGlobalElementId() {
        return this.maxMyGlobalElementId;
    }
    
    public int getMaxGlobalElementId() {
        return this.maxGlobalElementId;
    }
    
    public int[] getMyGlobalElementIds() {
        if (this.myGlobalElements != null) {
            return this.myGlobalElements;
        }
        
        int[] generatedGlobalElementIds = null;
        if (this.distributedLinearly) {
            generatedGlobalElementIds = new int[this.numMyGlobalElements];
            for (int i=0; i < this.numMyGlobalElements; i++) {
                generatedGlobalElementIds[i] = i + this.minMyGlobalElementId;
            }
        }
        
        return generatedGlobalElementIds;
    }
    
    public boolean isDistributedGlobally() {
        return this.distributedGlobally;
    }
    
    public boolean isDistributedUniformly() {
        return  this.distributedUniformly;
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
    
    protected void setComm(Comm comm) {
        this.comm = comm;
    }
    
    public void readExternal(ObjectInput in) throws java.io.IOException, ClassNotFoundException {
        this.numGlobalElements = in.readInt();
        this.numMyGlobalElements = in.readInt();
        this.minMyGlobalElementId = in.readInt();
        this.maxMyGlobalElementId = in.readInt();
        this.minLocalElementId = in.readInt();
        this.maxLocalElementId = in.readInt();
        this.minGlobalElementId = in.readInt();
        this.maxGlobalElementId = in.readInt();
        this.distributedGlobally = in.readBoolean();
        this.distributedLinearly = in.readBoolean();
        this.myGlobalElements = (int[]) in.readObject();
        this.gidsToLids = (TreeMap) in.readObject();
    }
    
    public void writeExternal(ObjectOutput out) throws java.io.IOException {
        out.writeInt(this.numGlobalElements);
        out.writeInt(this.numMyGlobalElements);
        out.writeInt(this.minMyGlobalElementId);
        out.writeInt(this.maxMyGlobalElementId);
        out.writeInt(this.minLocalElementId);
        out.writeInt(this.maxLocalElementId);
        out.writeInt(this.minGlobalElementId);
        out.writeInt(this.maxGlobalElementId);
        out.writeBoolean(this.distributedGlobally);
        out.writeBoolean(this.distributedLinearly);
        out.writeObject(this.myGlobalElements);
        out.writeObject(this.gidsToLids);
    }
    
    public boolean equals(ElementSpace otherElementSpace) {
        // see if they have the same reference address
        if (otherElementSpace == this) {
            return true;
        }
        
        // values that describe the distribution of Global and Local Elements
        if (otherElementSpace.numGlobalElements != this.numGlobalElements ||
            otherElementSpace.minGlobalElementId != this.minGlobalElementId ||
            otherElementSpace.maxGlobalElementId != this.maxGlobalElementId ||
            otherElementSpace.distributedGlobally != this.distributedGlobally ||
            otherElementSpace.distributedLinearly != this.distributedLinearly ||
            otherElementSpace.numMyGlobalElements != this.numMyGlobalElements ||
            otherElementSpace.minMyGlobalElementId != this.minMyGlobalElementId ||
            otherElementSpace.maxMyGlobalElementId != this.maxMyGlobalElementId ||
            otherElementSpace.minLocalElementId != this.minLocalElementId ||
            otherElementSpace.maxLocalElementId != this.maxLocalElementId) {
                
            /*this.println("STD", "numGlobalElements: " + otherElementSpace.numGlobalElements + " " + this.numGlobalElements);
            this.println("STD", "minGlobalElementId: " + otherElementSpace.minGlobalElementId + " " + this.minGlobalElementId);
            this.println("STD", "maxGlobalElementId: " + otherElementSpace.maxGlobalElementId + " " + this.maxGlobalElementId);
            this.println("STD", "distributedGlobally: " + otherElementSpace.distributedGlobally + " " + this.distributedGlobally);
            this.println("STD", "distributedLinearly: " + otherElementSpace.distributedLinearly + " " + this.distributedLinearly);
            this.println("STD", "numMyGlobalElements: " + otherElementSpace.numMyGlobalElements + " " + this.numMyGlobalElements);
            this.println("STD", "minMyGlobalElementId: " + otherElementSpace.minMyGlobalElementId + " " + this.minMyGlobalElementId);
            this.println("STD", "maxMyGlobalElementId: " + otherElementSpace.maxMyGlobalElementId + " " + this.maxMyGlobalElementId);
            this.println("STD", "minLocalElementId: " + otherElementSpace.minLocalElementId + " " + this.minLocalElementId);
            this.println("STD", "maxLocalElementId: " + otherElementSpace.maxLocalElementId + " " + this.maxLocalElementId);
            this.println("STD", "In ElementSpace.equals: one of the attributes of the two ElementSpaces did not match.");*/
            return false;
        }
        
        if (!this.distributedLinearly && !this.distributedUniformly && !otherElementSpace.gidsToLids.equals(this.gidsToLids)) {
            //this.println("STD", "In ElementSpace.equals: !otherElementSpace.gidsToLids.equals(this.gidsToLids)");
            return false;
        }
        
        return true;
    }
}
