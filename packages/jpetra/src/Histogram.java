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

package JpetraVis;

/* Histogram Printing Program
   Created by Kara Hansen
   Last modified: 3/19/02
*/

import javax.swing.*;
import java.awt.*;
import Jpetra.*;

public class Histogram{
    public static void main (String args[]) {

        //double vector[] = {0, 403, 487, -73, -87, -42, 34, 125, 264, 367, 320, 560, 104, 69, 0, 75, 98, 197, 300, 173, 982, 220, 200, 187, 95, 1075, 874, 736, 85, 80, 75, 70, 65, 64, 63, 62, 61, 1249, 1132, 984, 2304, 1500, 2874, 2385, 873, 762, 843, 127, 982, 763, 752, 99, 874, 370, 177, 199, 200, 564, 872, 1065, 1136, 1287, 1299, 673, 934, 185, 982, 543, 184, 67, 87, 107, 115, 60, 90, 70, 140, 80, 55, 43, 123, 345, 66, 76, 77, 143, 126, 55, 260, 197, 121, 96, 104, 52, 176, 199, 99};
        //double vector[] = {-4, 0, 3, 15, 23, 15, 1, 2, 3, 10, -2, -1, 5, 7, 24, 10, 22, 6, -1, -9, 0, 40, 33, 30, 25, 20, 34, 22, 19, 18, 13, 15, 24};
	double vector[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 9, 8, 5, 3, 2, 0, 4, 15, 12, 3, 4, 5, 9};
       
	VecView myVectorView = new VecView(vector, "Test");
	myVectorView.histogram();
	myVectorView.xygraph();    
    }
}
