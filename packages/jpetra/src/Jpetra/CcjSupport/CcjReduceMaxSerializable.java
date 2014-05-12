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

package Jpetra.CcjSupport;

import CCJ.*;
import java.io.Serializable;

/**
 * This class defines the algorithm used by by CCJ to compare all the <code>partialMins</code> during a CCJ
 * <code>allReduce</code> call.
 * Implements the CCJ interface <code>Reducible</code>.
 *
 * @author Jason Cross
 */

class CcjReduceMaxSerializable implements Reducible {

    public Serializable reduce(Serializable partialMax1, Serializable partialMax2) {
        if (partialMax1 == null) {
            return partialMax2;
        }
        if (partialMax2 == null) {
            return partialMax1;
        }
        
        Comparable partialMaxA = (Comparable) partialMax1;
        Comparable partialMaxB = (Comparable) partialMax2;
        
        // modifying paritalMaxA will modify the object passed in by the user
        // therefore for saftey reasons do not modify partialMaxA
            if (partialMaxB.compareTo(partialMaxA) < 0) {
                partialMaxB = partialMaxA;
        }
        
        return partialMax2;
    }
}
