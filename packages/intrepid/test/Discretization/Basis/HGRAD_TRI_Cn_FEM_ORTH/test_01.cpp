// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert Kirby (robert.c.kirby@ttu.edu)
//
// ************************************************************************
// @HEADER


/** \file
\brief  Unit test of Dubiner basis class
\author Created by R. Kirby
*/

#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Intrepid_HGRAD_TRI_Cn_FEM_ORTH.hpp"
#include "Intrepid_CubatureDirectTriDefault.hpp"
#include "Intrepid_PointTools.hpp"
#include "Shards_CellTopology.hpp"
#include <iostream>
using namespace Intrepid;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);
  
  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);
  
  *outStream \
    << "===============================================================================\n" \
    << "|                                                                             |\n" \
    << "|                           Unit Test OrthogonalBases                         |\n" \
    << "|                                                                             |\n" \
    << "|     1) Tests orthogonality of triangular orthogonal basis (Dubiner)         |\n" \
    << "|                                                                             |\n" \
    << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
    << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
    << "|                      Robert Kirby (robert.c.kirby@ttu.edu)                  |\n" \
    << "|                                                                             |\n" \
    << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
    << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
    << "|                                                                             |\n" \
    << "===============================================================================\n";
  
  int errorFlag  = 0;
  
  // First, get a reference quadrature rule

  CubatureDirectTriDefault<double,FieldContainer<double> > myCub(20);
  FieldContainer<double> cubPts( myCub.getNumPoints() , 2 );
  FieldContainer<double> cubWts( myCub.getNumPoints() );

  myCub.getCubature( cubPts , cubWts );
  
  // Tabulate the basis functions at the cubature points
  const int deg =3;
  const int polydim = (deg+1)*(deg+2)/2;
  FieldContainer<double> basisAtCubPts( polydim , myCub.getNumPoints() );
  
  Basis_HGRAD_TRI_Cn_FEM_ORTH<double,FieldContainer<double> > myBasis( deg );

  myBasis.getValues( basisAtCubPts , cubPts , OPERATOR_VALUE );

  // Now let's compute the mass matrix
  for (int i=0;i<polydim;i++) {
    for (int j=i;j<polydim;j++) {
      double cur = 0;
      for (int k=0;k<myCub.getNumPoints();k++) {
	cur += cubWts(k) * basisAtCubPts( i , k ) * basisAtCubPts( j , k );
      }
      if (i != j && fabs( cur ) > INTREPID_TOL) {
	errorFlag++;
      }
      else if (i == j && fabs( cur ) < INTREPID_TOL ) {
	errorFlag++;
      }

    }
  }

  // compare the points against FIAT-tabulated values on a lattice 
  shards::CellTopology myTri_3( shards::getCellTopologyData< shards::Triangle<3> >() );  
  const int np_lattice = PointTools::getLatticeSize( myTri_3 , deg , 0 );
  FieldContainer<double> lattice( np_lattice , 2);
  PointTools::getLattice<double,FieldContainer<double> >( lattice , 
							  myTri_3 , 
							  deg , 
							  0 , 
							  POINTTYPE_EQUISPACED );	 
				 
  FieldContainer<double> dBasisAtLattice( polydim , np_lattice , 2 );
  myBasis.getValues( dBasisAtLattice , lattice , OPERATOR_D1 );

  double fiat_vals[] =
    {
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      2.000000000000000e+00,
      1.000000000000000e+00,
      2.000000000000000e+00,
      1.000000000000000e+00,
      2.000000000000000e+00,
      1.000000000000000e+00,
      2.000000000000000e+00,
      1.000000000000000e+00,
      2.000000000000000e+00,
      1.000000000000000e+00,
      2.000000000000000e+00,
      1.000000000000000e+00,
      2.000000000000000e+00,
      1.000000000000000e+00,
      2.000000000000000e+00,
      1.000000000000000e+00,
      2.000000000000000e+00,
      1.000000000000000e+00,
      2.000000000000000e+00,
      1.000000000000000e+00,
      0.000000000000000e+00,
      3.000000000000000e+00,
      0.000000000000000e+00,
      3.000000000000000e+00,
      0.000000000000000e+00,
      3.000000000000000e+00,
      0.000000000000000e+00,
      3.000000000000000e+00,
      0.000000000000000e+00,
      3.000000000000000e+00,
      0.000000000000000e+00,
      3.000000000000000e+00,
      0.000000000000000e+00,
      3.000000000000000e+00,
      0.000000000000000e+00,
      3.000000000000000e+00,
      0.000000000000000e+00,
      3.000000000000000e+00,
      0.000000000000000e+00,
      3.000000000000000e+00,
      -6.000000000000000e+00,
      -2.000000000000000e+00,
      -2.000000000000000e+00,
      0.000000000000000e+00,
      2.000000000000000e+00,
      2.000000000000000e+00,
      6.000000000000000e+00,
      4.000000000000000e+00,
      -4.000000000000000e+00,
      -1.333333333333333e+00,
      -3.330669073875470e-16,
      6.666666666666665e-01,
      3.999999999999999e+00,
      2.666666666666666e+00,
      -2.000000000000000e+00,
      -6.666666666666667e-01,
      2.000000000000000e+00,
      1.333333333333333e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      -2.000000000000000e+00,
      -6.000000000000000e+00,
      -2.000000000000000e+00,
      -2.666666666666667e+00,
      -2.000000000000000e+00,
      6.666666666666663e-01,
      -2.000000000000000e+00,
      4.000000000000000e+00,
      1.333333333333333e+00,
      -2.666666666666667e+00,
      1.333333333333333e+00,
      6.666666666666663e-01,
      1.333333333333333e+00,
      3.999999999999999e+00,
      4.666666666666666e+00,
      6.666666666666661e-01,
      4.666666666666666e+00,
      3.999999999999999e+00,
      8.000000000000000e+00,
      4.000000000000000e+00,
      -0.000000000000000e+00,
      -8.000000000000000e+00,
      -0.000000000000000e+00,
      -8.000000000000000e+00,
      -0.000000000000000e+00,
      -7.999999999999999e+00,
      -0.000000000000000e+00,
      -8.000000000000000e+00,
      0.000000000000000e+00,
      -1.333333333333334e+00,
      0.000000000000000e+00,
      -1.333333333333334e+00,
      0.000000000000000e+00,
      -1.333333333333334e+00,
      0.000000000000000e+00,
      5.333333333333332e+00,
      0.000000000000000e+00,
      5.333333333333332e+00,
      0.000000000000000e+00,
      1.200000000000000e+01,
      1.200000000000000e+01,
      3.000000000000000e+00,
      -1.333333333333333e+00,
      -1.666666666666667e+00,
      -1.333333333333334e+00,
      3.333333333333327e-01,
      1.200000000000000e+01,
      9.000000000000000e+00,
      5.333333333333335e+00,
      1.333333333333334e+00,
      -1.333333333333333e+00,
      -6.666666666666670e-01,
      5.333333333333331e+00,
      3.999999999999999e+00,
      1.333333333333334e+00,
      3.333333333333335e-01,
      1.333333333333333e+00,
      9.999999999999997e-01,
      0.000000000000000e+00,
      0.000000000000000e+00,
      6.000000000000000e+00,
      9.000000000000000e+00,
      2.000000000000000e+00,
      -2.333333333333333e+00,
      -2.000000000000000e+00,
      -4.333333333333333e+00,
      -6.000000000000000e+00,
      3.000000000000000e+00,
      -5.333333333333333e+00,
      1.333333333333333e+00,
      -4.440892098500626e-16,
      -6.666666666666672e-01,
      5.333333333333332e+00,
      6.666666666666663e+00,
      -7.333333333333333e+00,
      -1.666666666666667e+00,
      7.333333333333331e+00,
      5.666666666666665e+00,
      0.000000000000000e+00,
      0.000000000000000e+00,
      2.000000000000000e+00,
      1.300000000000000e+01,
      2.000000000000000e+00,
      5.000000000000002e+00,
      2.000000000000000e+00,
      -2.999999999999999e+00,
      2.000000000000000e+00,
      -1.100000000000000e+01,
      -1.333333333333333e+00,
      -2.000000000000000e+00,
      -1.333333333333333e+00,
      -6.666666666666669e-01,
      -1.333333333333333e+00,
      6.666666666666661e-01,
      4.666666666666665e+00,
      -3.000000000000001e+00,
      4.666666666666665e+00,
      7.666666666666663e+00,
      2.000000000000000e+01,
      1.000000000000000e+01,
      0.000000000000000e+00,
      1.500000000000000e+01,
      0.000000000000000e+00,
      1.500000000000000e+01,
      0.000000000000000e+00,
      1.500000000000000e+01,
      0.000000000000000e+00,
      1.500000000000000e+01,
      -0.000000000000000e+00,
      -3.333333333333333e+00,
      -0.000000000000000e+00,
      -3.333333333333333e+00,
      -0.000000000000000e+00,
      -3.333333333333333e+00,
      0.000000000000000e+00,
      1.666666666666666e+00,
      0.000000000000000e+00,
      1.666666666666666e+00,
      0.000000000000000e+00,
      3.000000000000000e+01
    };
  
  int fiat_index_cur = 0;
  for (int i=0;i<polydim;i++) {
    for (int j=0;j<np_lattice;j++) {
      for (int k=0;k<2;k++) {
	if (std::abs( dBasisAtLattice(i,j,k) - fiat_vals[fiat_index_cur] ) > INTREPID_TOL ) {
	  errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          
          // Output the multi-index of the value where the error is:
          *outStream << " At multi-index { ";
          *outStream << i << " " << j << " " << k;
          *outStream << "}  computed value: " << dBasisAtLattice(i,j,k)
                     << " but correct value: " << fiat_vals[fiat_index_cur] << "\n";
	  *outStream << "Difference: " << std::abs( dBasisAtLattice(i,j,k) - fiat_vals[fiat_index_cur] ) << "\n";
        }
	fiat_index_cur++;
      }
    }
  }

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";
  
  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);
  
  return errorFlag;
}
