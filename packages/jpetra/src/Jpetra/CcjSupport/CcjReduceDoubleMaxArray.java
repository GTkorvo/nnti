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

package Jpetra.CcjSupport;

import CCJ.*;
import java.io.Serializable;

/**
 * This class defines the algorithm used by by CCJ to compare all the <code>partialMaxs</code> during a CCJ
 * <code>allReduce</code> call.
 * Implements the CCJ interface <code>Reducible</code>.
 *
 * @author Jason Cross
 */

class CcjReduceDoubleMaxArray implements Reducible {

    public Serializable reduce(Serializable partialMaxs1, Serializable partialMaxs2) {
        double [] partialMaxsA = (double []) partialMaxs1;
        double [] partialMaxsB = (double []) partialMaxs2;
    
        for(int i=0; i < partialMaxsA.length; i++) {
            // modifying paritalMaxsA will modify the array passed in by the user
            // therefore for saftey reasons do not modify partialMaxsA
            if (partialMaxsB[i] < partialMaxsA[i]) {
                partialMaxsB[i] = partialMaxsA[i];
            }
        }
        
        return partialMaxsB;
    }
}
