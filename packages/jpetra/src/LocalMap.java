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
 * LocalMap.java
 *
 * Created on July 23, 2001, 2:14 PM
 */

package Jpetra;

/**
 *
 * @author  Michael William Boldt
 * @version 
 */
public class LocalMap extends Jpetra.Map {

    /** Creates new LocalMap */
    public LocalMap(int numMyElements, int indexBase, Comm comm) {
        super(numMyElements, numMyElements, indexBase, comm);
        if(checkInput() != 0) {
            System.out.println("Replicated Local Map not the same size onf all PEs");
            System.exit(1);
        }
    }

    public LocalMap(LocalMap map) {
        super(map);
        isDistributedGlobal = false;
        if(checkInput() != 0) {
            System.out.println("Replicated Local Map not the same size on all PEs");
            System.exit(1);
        }
    }
    
    private int checkInput() {
        isDistributedGlobal = false;
        int [] tmp = new int [2];
        int [] res = new int [2];
        tmp[0] = numMyElements;
        tmp[1] = - numMyElements;
        res = comm.maxAll(tmp);
    
        int tmp1 = res[0];
        int tmp2 = - res[1];
    
        if(tmp1 == tmp2) return 0;
        else return -1;
    }
}
