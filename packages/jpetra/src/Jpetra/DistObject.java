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
 * All objects which can be constructed and distributed globally among the vnodes should extend this class.  <code>DistObject</code> presents a common interface to the user to use for all objects that can be distributed globally.  Thus, some methods must be overriden by classes extending <code>DistObject</code> while others must not be overridden.  <code>DistObject</code> also defines constants that correspond to the different avaliable combine modes.
 *
 * @author Jason Cross
 */
public abstract class DistObject extends JpetraObject {
    /**
     * This combine mode will cause all imported values to be added to existing values.
     */    
    public static final int ADD = 0;
    /**
     * This combine mode will only copy and permute local values.
     */    
    public static final int ZERO = 1;
    /**
     * This combine mode will cause all imported values to be averaged with existing values.
     */    
    public static final int AVERAGE = 2;
    /**
     * This combine mode replaces values with the maximum absolute value of either the imported value or the existing corresponding local value.
     */    
    public static final int ABSMAX = 3;
    
    private Distributor distributor;
    
    // for data I will receive
    private int numSameGids;
    private int[] remoteLids;  // another vnode's lid for the gid I need
    private int[] remoteGids;  // gid I need not on my vnode
    private int[] permuteToLids;
    private int[] permuteFromLids;
    private int combineMode;
    private DistObject distObjectSource;
    
    // for data I will send
    private int[] exportLids;
    private int[] exportGids;
    private int[] exportVnodeIds;
    
    /**
     * Does nothing.
     */    
    public DistObject() {
    }
    
    /**
     * Imports values from the source <code>DistObject</code> to this <code>DistObject<code> using the specified importer and combine mode.
     * This method should not be overriden by classes extending <code>DistObject</code>.
     *
     * @param distObjectSource The <code>DistObject</code> to import values from.
     * @param importer The <code>Import</code> object built from the source <code>VectorSpace</code> and the target <code>VectorSpace</code>.
     * @param combineMode One of the declared combine mode constants in <code>DistObject</code>.
     */    
    public void importValues(DistObject distObjectSource, Import importer, int combineMode) {
        this.numSameGids = importer.getNumSameGids();
        this.remoteLids = importer.getRemoteLids();
        this.remoteGids = importer.getRemoteGids();
        this.permuteToLids = importer.getPermuteToLids();
        this.permuteFromLids = importer.getPermuteFromLids();
        this.exportLids = importer.getExportLids();
        this.exportGids = importer.getExportGids();
        this.exportVnodeIds = importer.getExportVnodeIds();
        this.distributor = importer.getDistributor();
        
        this.combineMode = combineMode;
        this.distObjectSource = distObjectSource;
        
        doTransfer();
    }
    
    public void importValues(DistObject distObjectSource, Export exporter, int combineMode) {
        
    }
    
    public void exportValues(DistObject distObjectSource, Import importer, int combineMode) {
        
    }
    
    public void exportValues(DistObject distObjectSource, Export exporter, Comm comm, int combineMode) {
        
    }
    
    /**
     * Should not be called directly by the user.
     *
     * Calls <code>copyAndPermute</code>, <code>packAndPrepare</code>, and <code>unpackAndCombine</code> which are all
     * overridden by the calling child class.  <code>doTransfer</code> also calls <code>Distributor.distribute</code> in
     * order to send all exports to the corresponding vnodes and receive all imports from other vnodes.
     * This method should not be overriden by classes extending <code>DistObject</code>.
     */    
    public void doTransfer() {
        copyAndPermute(this.distObjectSource, this.numSameGids, this.permuteToLids, this.permuteFromLids, this.combineMode);
        if (this.combineMode == DistObject.ZERO) {
            // just doing a local copy and permute so we're done
            return;
        }
        Serializable[] exportData = packAndPrepare(this.distObjectSource, this.exportGids, this.exportLids);
        Serializable[] importData = this.distributor.distribute(exportData);
        unpackAndCombine(importData, combineMode);
    }
    
    /**
     * Should not be called directly by the user.
     *
     * This method must be overridden by classes extending <code>DistObject</code>.
     *
     * @param distObjectSource The object from which elements will be packedup for exportation to other vnodes.
     * @param exportGids The global IDs of the elements to be exported.
     * @param exportLids The local IDs of elements to be exported.
     *
     * @return The packed data that <code>Distributor.distribute</code> will export to other vnodes.
     */    
    public Serializable[] packAndPrepare(DistObject distObjectSource, int[] exportGids, int[] exportLids) {
        JpetraObject.println("ERR", "A fatal error has occurred.  packAndPrepare was not overridden, please contact the developer.");
        System.exit(1);
        return null;
    }
    
    /**
     * Should not be called directly by the user.
     *
     * This method must be overridden by classes extending <code>DistObject</code>.
     *
     * @param importData The data object returned from <code>Distributor.distribute</code> which contains all elements sent to this vnode.
     * @param combineMode One of the declared combine mode constants in <code>DistObject</code>.
     */    
    public void unpackAndCombine(Serializable[] importData, int combineMode) {
        this.println("ERR", "A fatal error has occurred.  unpackAndCombine was not overridden, please contact the developer.");
        System.exit(1);
    }
    
    /**
     * Should not be called directly by the user.
     *
     * This method must be overridden by classes extending <code>DistObject</code>.
     *
     * @param distObjectSource The object from which local values will be copied and permuted.
     * @param numSameGids The number of global ids that are the same starting from myMinGlobalEntryId.
     * @param permuteToLids The local IDs of this <code>DistObject</code> that will permute elements from <code>distObjectSource</code> according to the corresponding local IDs in <code>permuteFromLids</code>.
     * @param permuteFromLids The local IDs from <code>distObjectSource</code> that will be permuted to the local IDs in <code>permuteToLids</code> correspndong to this <code>DistObject</code>.
     * @param combineMode One of the declared combine mode constants in <code>DistObject</code>.
     */    
    public void copyAndPermute(DistObject distObjectSource, int numSameGids, int[] permuteToLids, int[] permuteFromLids, int combineMode) {
        this.println("ERR", "A fatal error has occurred.  copyAndPermute was not overridden, please contact the developer.");
        System.exit(1);
    }
}
