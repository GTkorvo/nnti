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

/*
 * CcjGatherDoubleArray.java
 *
 * Created on June 20, 2003, 7:58 AM
 */

package Jpetra.CcjSupport;

import CCJ.*;
import java.io.Serializable;

/**
 * This class defines the algorithm used by by CCJ to sum all the <code>partialSums</code> together during a CCJ
 * <code>allReduce</code> call.
 * Implements the CCJ interface <code>Reducible</code>.
 *
 * @author Jason Cross
 */

class CcjReduceAllSumIntArray implements Reducible {

    public Serializable reduce(Serializable partialSums1, Serializable partialSums2) {

        //cast the partialSums into a useful object type 
        int [] partialSumsA = (int []) partialSums1;
        int [] partialSumsB = (int []) partialSums2;
    
        //add the partial sums together
        for(int i=0; i < partialSumsA.length; i++) {
            partialSumsB[i] += partialSumsA[i];
        }
        
        //returned the summed partial sums
        return partialSumsB;
    }
}
