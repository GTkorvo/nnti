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
 * Objects that extend <code>DistObject</code> can be constructed and distributed globally among the vnodes.  <code>DistObject</code> presents a common interface to the user to use for all objects that can be distributed globally.  Thus, some methods must be overriden by classes extending <code>DistObject</code> while others must not be overridden.  <code>DistObject</code> also defines constants that correspond to the different avaliable combine modes.
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
    /**
     * This combine mode replaces existing values with the imported ones.
     */
    public static final int REPLACE = 4;
    
    private Distributor distributor;
    
    // for data I will receive
    //private int numSameGids;
    //private int[] remoteLids;  // another vnode's lid for the gid I need
    //private int[] remoteGids;  // gid I need not on my vnode
    //private int[] permuteToLids;
    //private int[] permuteFromLids;
    //private int combineMode;
    //private DistObject distObjectSource;
    
    // for data I will send
    //private int[] exportLids;
    //private int[] exportGids;
    //private int[] exportVnodeIds;
    
    // for flops
    private FlopCounter flops;
    
    /**
     * Does nothing.
     */
    public DistObject() {
        // empty
    }
    
    /**
     * Imports values from the source <code>DistObject</code> to this <code>DistObject</code> using the specified importer and combine mode.
     * This method should not be overriden by classes extending <code>DistObject</code>.
     *
     * @param distObjectSource The <code>DistObject</code> to import values from.
     * @param importer The <code>Import</code> object built from the source <code>VectorSpace</code> and the target <code>VectorSpace</code>.
     * @param combineMode One of the declared combine mode constants in <code>DistObject</code>.
     */
    public void importValues(DistObject distObjectSource, Import importer, int combineMode) {
        this.distributor = importer.getDistributor();
        
        
        /*
        this.remoteLids = importer.getRemoteLids();
        this.remoteGids = importer.getRemoteGids();
        this.exportVnodeIds = importer.getExportVnodeIds();*/
        
        //this.numSameGids = importer.getNumSameGids();
        //this.permuteToLids = importer.getPermuteToLids();
        //this.permuteFromLids = importer.getPermuteFromLids();
        //this.exportLids = importer.getExportLids();
        //this.exportGids = importer.getExportGids();
        
        
        
        //this.combineMode = combineMode;
        //this.distObjectSource = distObjectSource;
        
        doTransfer(distObjectSource, importer.getNumSameGids(), importer.getPermuteToLids(), importer.getPermuteFromLids(),
                   importer.getExportGids(), importer.getExportLids(), combineMode, false);
    }
    
    public void importValues(DistObject distObjectSource, Export exporter, int combineMode) {

        this.distributor = exporter.getDistributor();
        if ((!distObjectSource.getVectorSpace().getComm().isSerial()) && (this.distributor.doneForwardOp() == false)) {
            // do a null forward op to figure out what gids will be sent to us
            // so we know which gids we must send for the reverse op
            NullDistObject nullDistObjectTarget = new NullDistObject(this.getVectorSpace());
            NullDistObject nullDistObjectSource = new NullDistObject(distObjectSource.getVectorSpace());
            //this.println("STD", "Doing null export...");
            nullDistObjectSource.exportValues(nullDistObjectTarget, exporter, DistObject.ADD);
            //this.println("STD", "Done doing null export...");
            // the null forward op is complete, now we can do the reverse op
        }
        
        // same as forward op
        //this.numSameGids = exporter.getNumSameGids();
        
        // reverse of forward op
        //this.remoteLids = exporter.getExportLids();
        //this.remoteGids = exporter.getExportGids();
        //this.permuteToLids = exporter.getPermuteFromLids();
        //this.permuteFromLids = exporter.getPermuteToLids();
        if (!distObjectSource.getVectorSpace().getComm().isSerial()) {
        //    this.exportLids = this.distributor.getReverseExportLids();
        //    this.exportGids = this.distributor.getReverseExportGids();
        //    this.exportVnodeIds = this.distributor.getReverseExportVnodeIds();
        }
        
        //this.combineMode = combineMode;
        //this.distObjectSource = distObjectSource;
        
       doTransfer(distObjectSource, exporter.getNumSameGids(), exporter.getPermuteToLids(), exporter.getPermuteFromLids(),
                   exporter.getExportGids(), exporter.getExportLids(), combineMode, true);
    }
    
    public void exportValues(DistObject distObjectSource, Import importer, int combineMode) {
        this.distributor = importer.getDistributor();
        if ((!distObjectSource.getVectorSpace().getComm().isSerial()) && (this.distributor.doneForwardOp() == false)) {
            // do a null forward op to figure out what gids will be sent to us
            // so we know which gids we must send for the reverse op
            NullDistObject nullDistObjectTarget = new NullDistObject(this.getVectorSpace());
            NullDistObject nullDistObjectSource = new NullDistObject(distObjectSource.getVectorSpace());
            //this.println("STD", "Doing null export...");
            nullDistObjectSource.importValues(nullDistObjectTarget, importer, DistObject.ADD);
            //this.println("STD", "Done doing null export...");
            // the null forward op is complete, now we can do the reverse op
        }
        
        // same as forward op
        //this.numSameGids = importer.getNumSameGids();
        
        // reverse of forward op
        //this.remoteLids = importer.getExportLids();
        //this.remoteGids = importer.getExportGids();
        //this.permuteToLids = importer.getPermuteFromLids();
        //this.permuteFromLids = importer.getPermuteToLids();
        if (!distObjectSource.getVectorSpace().getComm().isSerial()) {
          //  this.exportLids = this.distributor.getReverseExportLids();
          //  this.exportGids = this.distributor.getReverseExportGids();
          //  this.exportVnodeIds = this.distributor.getReverseExportVnodeIds();
        }
        
        //this.combineMode = combineMode;
        //this.distObjectSource = distObjectSource;
        doTransfer(distObjectSource, importer.getNumSameGids(), importer.getPermuteToLids(), importer.getPermuteFromLids(),
                   importer.getExportGids(), importer.getExportLids(), combineMode, true);
    }
    
    public void exportValues(DistObject distObjectSource, Export exporter, int combineMode) {
        
        //this.numSameGids = exporter.getNumSameGids();
        //this.remoteLids = exporter.getRemoteLids();
        //this.remoteGids = exporter.getRemoteGids();
        //this.permuteToLids = exporter.getPermuteToLids();
        //this.permuteFromLids = exporter.getPermuteFromLids();
        //this.exportLids = exporter.getExportLids();
        //this.exportGids = exporter.getExportGids();
        //this.exportVnodeIds = exporter.getExportVnodeIds();
        this.distributor = exporter.getDistributor();
        
        //this.combineMode = combineMode;
        //this.distObjectSource = distObjectSource;
        
        doTransfer(distObjectSource, exporter.getNumSameGids(), exporter.getPermuteToLids(), exporter.getPermuteFromLids(),
                   exporter.getExportGids(), exporter.getExportLids(), combineMode, false);
    }
    
    /**
     * Should not be called directly by the user.
     *
     * Calls <code>copyAndPermute</code>, <code>packAndPrepare</code>, and <code>unpackAndCombine</code> which are all
     * overridden by the calling child class.  <code>doTransfer</code> also calls <code>Distributor.distribute</code> in
     * order to send all exports to the corresponding vnodes and receive all imports from other vnodes.
     * This method should not be overriden by classes extending <code>DistObject</code>.
     */
    public void doTransfer(DistObject distObjectSource, int numSameGids, int[] permuteToLids, int[] permuteFromLids, int[] exportGids, int[] exportLids, int combineMode, boolean doReverse) {
        copyAndPermute(distObjectSource, numSameGids, permuteToLids, permuteFromLids, combineMode);
        if (distObjectSource.getVectorSpace().getComm().isSerial() || (combineMode == DistObject.ZERO)) {
            // just doing a local copy and permute so we're done
            //this.println("STD", "DistObject: Doing copyAndPermute Only!");
            return;
        }
        Serializable[] exportData = packAndPrepare(distObjectSource, exportGids, exportLids);
        Serializable[] importData = this.distributor.distribute(exportData, doReverse);
        int[][] reverseExportVnodeIdsGidsLids = unpackAndCombine(importData, combineMode);
        if ((doReverse != true) && (this.distributor.doneForwardOp() != true)) {
            this.distributor.setDoneForwardOp(true);
            this.distributor.setReverseExportVnodeIdsGidsLids(reverseExportVnodeIdsGidsLids);
        }
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
    public abstract Serializable[] packAndPrepare(DistObject distObjectSource, int[] exportGids, int[] exportLids);
    
    /**
     * Should not be called directly by the user.
     *
     * This method must be overridden by classes extending <code>DistObject</code>.
     *
     * @param importData The data object returned from <code>Distributor.distribute</code> which contains all elements sent to this vnode.
     * @param combineMode One of the declared combine mode constants in <code>DistObject</code>.
     */
    public abstract int[][] unpackAndCombine(Serializable[] importData, int combineMode);
    
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
    public abstract void copyAndPermute(DistObject distObjectSource, int numSameGids, int[] permuteToLids, int[] permuteFromLids, int combineMode);

    public abstract VectorSpace getVectorSpace();
    
    public void setFlopCounter(FlopCounter flops) {
        this.flops = flops;
    }
    
    public void updateFlops(double numFlops) {
        if (this.flops != null) {
            this.flops.updateFlops(numFlops);
        }
    }
    
    public void resetFlops() {
        if (this.flops != null) {
            this.flops.resetFlops();
        }
    }
    
    public double getFlops() {
        if (this.flops != null) {
            return this.flops.getFlops();
        }
        
        return 0;
    }
    
    public double getGlobalMegaFlops(Comm comm) {
        if (this.flops != null) {
            double[] flopsDone = comm.sumAll(new double[]{this.flops.getFlops()});
            return flopsDone[0]/(1000000);
        }
        
        return 0;
    }
    
    public FlopCounter getFlopCounter() {
        return this.flops;
    }
    
}
