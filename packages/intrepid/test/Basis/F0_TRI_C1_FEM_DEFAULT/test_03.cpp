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
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file test_02.cpp
\brief  Unit tests for the Intrepid::LocalForm0  class,
        with Basis_F0_TRI_C1_FEM_DEFAULT basis.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_LocalForm0.hpp"
#include "Intrepid_DefaultBasisFactory.hpp"
#include "Intrepid_CubatureDirect.hpp"
#include "Intrepid_MultiCell.hpp"
#include "Intrepid_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Time.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
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
  << "|                 Unit Test (Basis_F0_TRI_C1_FEM_DEFAULT)                     |\n" \
  << "|                                                                             |\n" \
  << "|     1) Basis creation, computation of basis function values                 |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 1: Basic performance measurements                                      |\n"\
  << "===============================================================================\n";

  int errorFlag = 0;

  try {
    DefaultBasisFactory<double> BFactory;
    Teuchos::RCP<Basis<double> > basis =  BFactory.create(
      FIELD_FORM_0, CELL_TRI, RECONSTRUCTION_SPACE_COMPLETE, 1, BASIS_FEM_DEFAULT, COORDINATES_CARTESIAN);

    Teuchos::Array<Teuchos::Array< Teuchos::RCP<Cubature<double> > > > allCubs;
    Teuchos::RCP<Cubature<double> > cellCub = Teuchos::rcp(new CubatureDirect<double>(CELL_TRI,2) );
    allCubs.resize(1); 
    allCubs[0].resize(1);
    allCubs[0][0] = cellCub;
    /*************** use if you need face cubatures
    Teuchos::RCP<Cubature<double> > faceCub = Teuchos::rcp(new CubatureDirect<double>(CELL_EDGE,2) );
    allCubs.resize(2); 
    allCubs[0].resize(1);
    allCubs[0][0] = cellCub;
    allCubs[1].resize(3);
    allCubs[1][0] = faceCub; 
    allCubs[1][1] = faceCub; 
    allCubs[1][2] = faceCub;
    ****************/

    string basedir = "./data";

    int nCells = 100;
    double triNodes[3*2*100];
    int chunkSize = 4;
    int nChunks = nCells/chunkSize;
    int skip = 3*2*chunkSize;

    double triChunkNodes[] = {
      // tri 0
      0.0, 0.0,
      1.0, 0.0,
      0.0, 1.0,
      // tri 1
      4.0, 5.0,
      -6.0, 2.0,
      4.0, -3.0,
      // tri 2
      -6.0, -3.0,
      9.0, 2.0,
      8.9, 2.1,
      // tri 3
      -6.0, -3.0,
      12.0, 3.0,
      2.9, 0.1
    };

    for (int i=0; i<nChunks; i++) {
      for (int j=0; j<skip; j++) {
        triNodes[i*skip+j] = triChunkNodes[j];
      }
    }

    MultiCell<double> largeMCell(nCells,          // number of cells (triangles) in the multicell instance
                                 CELL_TRI,        // generating cell type
                                 triNodes);       // array with interleaved node coordinates

    LexContainer<double> massMatrices;

    for (ECompEngine compEng = COMP_CPP; compEng < COMP_ENGINE_MAX; compEng++) {
      Teuchos::Time timer("Timer - One Large Multicell");
      timer.start();
      for (int cubDeg=2; cubDeg<=20; cubDeg++) {
        // set cubature
        cellCub = Teuchos::rcp(new CubatureDirect<double>(CELL_TRI,cubDeg) );
        allCubs[0][0] = cellCub;
        // create local form
        LocalForm0<double> form0(basis, allCubs, compEng);
        // compute mass matrices
        form0.getOperator(massMatrices, OPERATOR_VALUE, OPERATOR_VALUE, largeMCell);
      }
      timer.stop();
      *outStream << "\n" << timer.name() << "\n";
      *outStream << "Computational engine: " << ECompEngineToString(compEng) << "\n";
      *outStream << "Measured time in sec: " << timer.totalElapsedTime() << "\n";
    }

    MultiCell<double> smallMCell(chunkSize,       // number of cells (triangles) in the multicell instance
                                 CELL_TRI,        // generating cell type
                                 triChunkNodes);  // array with interleaved node coordinates

    for (ECompEngine compEng = COMP_CPP; compEng < COMP_ENGINE_MAX; compEng++) {
      Teuchos::Time timer("Timer - Many Small Multicells");
      timer.start();
      for (int ex=0; ex<nChunks; ex++) {
        for (int cubDeg=2; cubDeg<=20; cubDeg++) {
          // set cubature
          cellCub = Teuchos::rcp(new CubatureDirect<double>(CELL_TRI,cubDeg) );
          allCubs[0][0] = cellCub;
          // create local form
          LocalForm0<double> form0(basis, allCubs, compEng);
          // compute mass matrices
          form0.getOperator(massMatrices, OPERATOR_VALUE, OPERATOR_VALUE, smallMCell);
        }
      }
      timer.stop();
      *outStream << "\n" << timer.name() << "\n";
      *outStream << "Computational engine: " << ECompEngineToString(compEng) << "\n";
      *outStream << "Measured time in sec: " << timer.totalElapsedTime() << "\n";
    }
  }
  catch (std::logic_error err) {
      *outStream << err.what() << "\n\n";
          errorFlag = -999;
  };
            

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
