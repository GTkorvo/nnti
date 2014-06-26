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

import java.io.Externalizable;
import java.io.ObjectOutput;
import java.io.ObjectInput;

/**
 *
 * @author  Jason Cross
 */
public class VectorSpace extends JpetraObject implements Externalizable {
    ElementSpace elementSpace;
    Directory directory;
    
    public VectorSpace() {
        // empty
    }
    
    public VectorSpace(ElementSpace elementSpace) {
        this.elementSpace = elementSpace;
    }
    
    public boolean isCompatible(VectorSpace otherVectorSpace) {
        if (otherVectorSpace == this) {
            return true;
        }
        if (this.getNumMyGlobalEntries() == otherVectorSpace.getNumMyGlobalEntries()) {
            return true;
        }
        
        return false;
    }
    
    public int getNumGlobalEntries() {
        return this.elementSpace.getNumGlobalElements();
    }
    
    public int getNumMyGlobalEntries() {
        return this.elementSpace.getNumMyGlobalElements();
    }
    
    public int getMinLocalIndex() {
        return this.elementSpace.getMinLocalElementId();
    }
    
    public int getMyMaxGlobalIndex() {
        return this.elementSpace.getMyMaxGlobalElementId();
    }
    
    public int getMyMinGlobalIndex() {
        return this.elementSpace.getMyMinGlobalElementId();
    }
    
    public int getLocalIndex(int globalIndex) {
        return this.elementSpace.getLocalElementId(globalIndex);
    }
    
    public int getGlobalIndex(int localIndex) {
        return this.elementSpace.getGlobalElementId(localIndex);
    }
    
    public boolean isMyLocalIndex(int localIndex) {
        return this.elementSpace.isMyLocalElementId(localIndex);
    }
    
    public boolean isMyGlobalIndex(int globalIndex) {
        return this.elementSpace.isMyGlobalElementId(globalIndex);
    }
    
    public int[] getMyGlobalEntryIds() {
        return this.elementSpace.getMyGlobalElementIds();
    }
    
    public int[][] getRemoteVnodeIdList(int[] remoteGlobalIndicies) {
        if (this.directory == null) {
            this.directory = new BasicDirectory(this);
        }
        
        return directory.getDirectoryEntries(remoteGlobalIndicies);
    }
    
    public boolean isDistributedGlobally() {
        return this.elementSpace.isDistributedGlobally();
    }
    
    public boolean isDistributedLinearly() {
        return this.elementSpace.isDistributedLinearly();
    }
    
    public boolean isDistributedUniformly() {
        return this.elementSpace.isDistributedUniformly();
    }
    
    public int getMinGlobalEntryId() {
        return this.elementSpace.getMinGlobalElementId();
    }
    
    public int getMaxGlobalEntryId() {
        return this.elementSpace.getMaxGlobalElementId();
    }
    
    public Comm getComm() {
        return this.elementSpace.getComm();
    }
    
    public int getNumRemainderIndices() {
        return this.elementSpace.getNumRemainderIndices();
    }
    
    public int getNumIndicesPerVnode() {
        return this.elementSpace.getNumIndicesPerVnode();
    }
    
    protected void setComm(Comm comm) {
        this.elementSpace.setComm(comm);
    }
    
    // !! not yet implemented
    public boolean equals(VectorSpace otherVectorSpace) {
        // check to see if they have the same reference address
        if (this == otherVectorSpace) {
            return true;
        }
        
        return otherVectorSpace.elementSpace.equals(this.elementSpace);
    }
    
    public void readExternal(ObjectInput in) throws java.io.IOException, ClassNotFoundException {
        this.elementSpace = (ElementSpace) in.readObject();
    }
    
    public void writeExternal(ObjectOutput out) throws java.io.IOException {
        out.writeObject(this.elementSpace);
    }
    
}
