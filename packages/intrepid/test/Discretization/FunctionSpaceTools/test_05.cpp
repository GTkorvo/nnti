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

/** \file test_05.cpp
\brief  Unit test for the FunctionSpaceTools class, testing H-div,
        so similar to test_03.cpp; however, shards::Arrays are
        used as computational arrays and wrapped as FieldContainers.

\author Created by D. Ridzal, P. Bochev, and K. Peterson.
*/

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Shards_Array.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using namespace std;
using namespace Intrepid;

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Cell )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Cell )

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Field )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Field )

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Point )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Point )

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Dim )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Dim )

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( Node )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( Node )


typedef FieldContainer<double> FC;


int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

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
  << "|                      Unit Test (FunctionSpaceTools)                         |\n" \
  << "|                                                                             |\n" \
  << "|     1) basic operator transformations and integration in HDIV               |\n" \
  << "|                    via shards::Array wrappers                               |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n";


  int errorFlag = 0;

  typedef FunctionSpaceTools fst; 

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 1: correctness of math operations                                      |\n"\
  << "===============================================================================\n";

  outStream->precision(20);

  try {
      shards::CellTopology cellType = shards::getCellTopologyData< shards::Hexahedron<> >();    // cell type: hex

      /* Related to cubature. */
      DefaultCubatureFactory<double> cubFactory;                                                // create cubature factory
      int cubDegree = 20;                                                                       // cubature degree
      Teuchos::RCP<Cubature<double> > myCub = cubFactory.create(cellType, cubDegree);           // create default cubature
      int spaceDim = myCub->getDimension();                                                     // get spatial dimension 
      int numCubPoints = myCub->getNumPoints();                                                 // get number of cubature points

      /* Related to basis. */
      Basis_HDIV_HEX_I1_FEM<double, FieldContainer<double> > hexBasis;                          // create H-div basis on a hex
      int numFields = hexBasis.getCardinality();                                                // get basis cardinality
 
      /* Cell geometries and orientations. */
      int numCells    = 4;
      int numNodes    = 8;
      int numCellData = numCells*numNodes*spaceDim;
      int numSignData = numCells*numFields;

      double hexnodes[] = {
        // hex 0  -- affine
        -1.0, -1.0, -1.0,
        1.0, -1.0, -1.0,
        1.0, 1.0, -1.0,
        -1.0, 1.0, -1.0,
        -1.0, -1.0, 1.0,
        1.0, -1.0, 1.0,
        1.0, 1.0, 1.0,
        -1.0, 1.0, 1.0,
        // hex 1  -- affine
        -3.0, -3.0, 1.0,
        6.0, 3.0, 1.0,
        7.0, 8.0, 0.0,
        -2.0, 2.0, 0.0,
        -3.0, -3.0, 4.0,
        6.0, 3.0, 4.0,
        7.0, 8.0, 3.0,
        -2.0, 2.0, 3.0,
        // hex 2  -- affine
        -3.0, -3.0, 0.0,
        9.0, 3.0, 0.0,
        15.0, 6.1, 0.0,
        3.0, 0.1, 0.0,
        9.0, 3.0, 0.1,
        21.0, 9.0, 0.1,
        27.0, 12.1, 0.1,
        15.0, 6.1, 0.1,
        // hex 3  -- nonaffine
        -2.0, -2.0, 0.0,
        2.0, -1.0, 0.0,
        1.0, 6.0, 0.0,
        -1.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        1.0, 0.0, 1.0,
        1.0, 1.0, 1.0,
        0.0, 1.0, 1.0
      };

      short facesigns[] = {
        1, 1, 1, 1, 1, 1,
        1, -1, 1, -1, 1, -1,
        -1, -1, 1, 1, -1, 1,
        -1, -1, 1, 1, -1, -1
      };

      /* Computational arrays. */
      // First allocate one very large work space.
      Teuchos::Array<double> work_space(numCubPoints*spaceDim +
                                        numCubPoints +
                                        numCells*numNodes*spaceDim +
                                        numCells*numCubPoints*spaceDim*spaceDim +
                                        2*numCells*numCubPoints +
                                        numFields*numCubPoints +
                                        2*numCells*numFields*numCubPoints +
                                        numCells*numFields*numFields +
                                        numFields*numCubPoints*spaceDim +
                                        2*numCells*numFields*numCubPoints*spaceDim +
                                        numCells*numFields*numFields
                                       );

      int offset = 0;
      shards::Array<double,shards::NaturalOrder,Point,Dim> cub_points(&work_space[offset], numCubPoints, spaceDim);
      FC cub_points_FC(cub_points);
      offset += numCubPoints*spaceDim;
      shards::Array<double,shards::NaturalOrder,Point> cub_weights(&work_space[offset], numCubPoints);
      FC cub_weights_FC(cub_weights);
      offset += numCubPoints;
      shards::Array<double,shards::NaturalOrder,Cell,Node,Dim> cell_nodes(&work_space[offset], numCells, numNodes, spaceDim);
      FC cell_nodes_FC(cell_nodes);
      offset += numCells*numNodes*spaceDim;
      FieldContainer<short>  field_signs(numCells, numFields);
      shards::Array<double,shards::NaturalOrder,Cell,Point,Dim,Dim> jacobian(&work_space[offset], numCells, numCubPoints, spaceDim, spaceDim);
      FC jacobian_FC(jacobian);
      offset += numCells*numCubPoints*spaceDim*spaceDim;
      //shards::Array<double,shards::NaturalOrder> jacobian_inv(&work_space[offset]numCells, numCubPoints, spaceDim, spaceDim);
      shards::Array<double,shards::NaturalOrder,Cell,Point> jacobian_det(&work_space[offset], numCells, numCubPoints);
      FC jacobian_det_FC(jacobian_det);
      offset += numCells*numCubPoints;
      shards::Array<double,shards::NaturalOrder,Cell,Point> weighted_measure(&work_space[offset], numCells, numCubPoints);
      FC weighted_measure_FC(weighted_measure);
      offset += numCells*numCubPoints;

      shards::Array<double,shards::NaturalOrder,Field,Point> div_of_basis_at_cub_points(&work_space[offset], numFields, numCubPoints);
      FC div_of_basis_at_cub_points_FC(div_of_basis_at_cub_points);
      offset += numFields*numCubPoints;
      shards::Array<double,shards::NaturalOrder,Cell,Field,Point>
        transformed_div_of_basis_at_cub_points(&work_space[offset], numCells, numFields, numCubPoints);
      FC transformed_div_of_basis_at_cub_points_FC(transformed_div_of_basis_at_cub_points);
      offset += numCells*numFields*numCubPoints;
      shards::Array<double,shards::NaturalOrder,Cell,Field,Point>
        weighted_transformed_div_of_basis_at_cub_points(&work_space[offset], numCells, numFields, numCubPoints);
      FC weighted_transformed_div_of_basis_at_cub_points_FC(weighted_transformed_div_of_basis_at_cub_points);
      offset += numCells*numFields*numCubPoints;
      shards::Array<double,shards::NaturalOrder,Cell,Field,Field> stiffness_matrices(&work_space[offset], numCells, numFields, numFields);
      FC stiffness_matrices_FC(stiffness_matrices);
      offset += numCells*numFields*numFields;

      shards::Array<double,shards::NaturalOrder,Field,Point,Dim> value_of_basis_at_cub_points(&work_space[offset], numFields, numCubPoints, spaceDim);
      FC value_of_basis_at_cub_points_FC(value_of_basis_at_cub_points);
      offset += numFields*numCubPoints*spaceDim;
      shards::Array<double,shards::NaturalOrder,Cell,Field,Point,Dim>
        transformed_value_of_basis_at_cub_points(&work_space[offset], numCells, numFields, numCubPoints, spaceDim);
      FC transformed_value_of_basis_at_cub_points_FC(transformed_value_of_basis_at_cub_points);
      offset += numCells*numFields*numCubPoints*spaceDim;
      shards::Array<double,shards::NaturalOrder,Cell,Field,Point,Dim>
        weighted_transformed_value_of_basis_at_cub_points(&work_space[offset], numCells, numFields, numCubPoints, spaceDim);
      FC weighted_transformed_value_of_basis_at_cub_points_FC(weighted_transformed_value_of_basis_at_cub_points);
      offset += numCells*numFields*numCubPoints*spaceDim;
      shards::Array<double,shards::NaturalOrder,Cell,Field,Field> mass_matrices(&work_space[offset], numCells, numFields, numFields);
      FC mass_matrices_FC(mass_matrices);


      /******************* START COMPUTATION ***********************/

      // get cubature points and weights
      myCub->getCubature(cub_points_FC, cub_weights_FC);

      // fill cell vertex array
      cell_nodes_FC.setValues(hexnodes, numCellData);

      // set basis function signs, for each cell
      field_signs.setValues(facesigns, numSignData);

      // compute geometric cell information
      CellTools<double>::setJacobian(jacobian_FC, cub_points_FC, cell_nodes_FC, cellType);
      //CellTools<double>::setJacobianInv(jacobian_inv, jacobian);
      CellTools<double>::setJacobianDet(jacobian_det_FC, jacobian_FC);

      // compute weighted measure
      fst::computeCellMeasure<double>(weighted_measure_FC, jacobian_det_FC, cub_weights_FC);


      // Computing stiffness matrices:
      // tabulate divergences of basis functions at (reference) cubature points
      hexBasis.getValues(div_of_basis_at_cub_points_FC, cub_points_FC, OPERATOR_DIV);

      // transform divergences of basis functions
      fst::HDIVtransformDIV<double>(transformed_div_of_basis_at_cub_points_FC,
                                    jacobian_det_FC,
                                    div_of_basis_at_cub_points_FC);

      // multiply with weighted measure
      fst::multiplyMeasure<double>(weighted_transformed_div_of_basis_at_cub_points_FC,
                                   weighted_measure_FC,
                                   transformed_div_of_basis_at_cub_points_FC);

      // compute stiffness matrices
      fst::integrate<double>(stiffness_matrices_FC,
                             transformed_div_of_basis_at_cub_points_FC,
                             weighted_transformed_div_of_basis_at_cub_points_FC,
                             COMP_CPP);

      // apply field signs
      fst::applyLeftFieldSigns<double>(stiffness_matrices_FC, field_signs);
      fst::applyRightFieldSigns<double>(stiffness_matrices_FC, field_signs);


      // Computing mass matrices:
      // tabulate values of basis functions at (reference) cubature points
      hexBasis.getValues(value_of_basis_at_cub_points_FC, cub_points_FC, OPERATOR_VALUE);

      // transform values of basis functions
      fst::HDIVtransformVALUE<double>(transformed_value_of_basis_at_cub_points_FC,
                                      jacobian_FC,
                                      jacobian_det_FC,
                                      value_of_basis_at_cub_points_FC);

      // multiply with weighted measure
      fst::multiplyMeasure<double>(weighted_transformed_value_of_basis_at_cub_points_FC,
                                   weighted_measure_FC,
                                   transformed_value_of_basis_at_cub_points_FC);

      // compute mass matrices
      fst::integrate<double>(mass_matrices_FC,
                             transformed_value_of_basis_at_cub_points_FC,
                             weighted_transformed_value_of_basis_at_cub_points_FC,
                             COMP_CPP);

      // apply field signs
      fst::applyLeftFieldSigns<double>(mass_matrices_FC, field_signs);
      fst::applyRightFieldSigns<double>(mass_matrices_FC, field_signs);

      /*******************  STOP COMPUTATION ***********************/


      /******************* START COMPARISON ***********************/
      string basedir = "./testdata";
      for (int cell_id = 0; cell_id < numCells-1; cell_id++) {

        stringstream namestream;
        string filename;
        namestream <<  basedir << "/mass_HDIV_HEX_I1_FEM" << "_" << "0" << cell_id+1 << ".dat";
        namestream >> filename;

        ifstream massfile(&filename[0]);
        if (massfile.is_open()) {
          if (compareToAnalytic<double>(&mass_matrices(cell_id, 0, 0), massfile, 1e-10, 1) > 0)
            errorFlag++;
          massfile.close();
        }
        else {
          errorFlag = -1;
          std::cout << "End Result: TEST FAILED\n";
          return errorFlag;
        }

        namestream.clear();
        namestream << basedir << "/stiff_HDIV_HEX_I1_FEM" << "_" << "0" << cell_id+1 << ".dat";
        namestream >> filename;

        ifstream stifffile(&filename[0]);
        if (stifffile.is_open())
        {
          if (compareToAnalytic<double>(&stiffness_matrices(cell_id, 0, 0), stifffile, 1e-10, 1) > 0)
            errorFlag++;
          stifffile.close();
        }
        else {
          errorFlag = -1;
          std::cout << "End Result: TEST FAILED\n";
          return errorFlag;
        }

      }
      for (int cell_id = 3; cell_id < numCells; cell_id++) {

        stringstream namestream;
        string filename;
        namestream <<  basedir << "/mass_fp_HDIV_HEX_I1_FEM" << "_" << "0" << cell_id+1 << ".dat";
        namestream >> filename;

        ifstream massfile(&filename[0]);
        if (massfile.is_open()) {
          if (compareToAnalytic<double>(&mass_matrices(cell_id, 0, 0), massfile, 1e-4, 1, INTREPID_UTILS_SCALAR) > 0)
            errorFlag++;
          massfile.close();
        }
        else {
          errorFlag = -1;
          std::cout << "End Result: TEST FAILED\n";
          return errorFlag;
        }

        namestream.clear();
        namestream << basedir << "/stiff_fp_HDIV_HEX_I1_FEM" << "_" << "0" << cell_id+1 << ".dat";
        namestream >> filename;

        ifstream stifffile(&filename[0]);
        if (stifffile.is_open())
        {
          if (compareToAnalytic<double>(&stiffness_matrices(cell_id, 0, 0), stifffile, 1e-4, 1, INTREPID_UTILS_SCALAR) > 0)
            errorFlag++;
          stifffile.close();
        }
        else {
          errorFlag = -1;
          std::cout << "End Result: TEST FAILED\n";
          return errorFlag;
        }

      }
      /******************* STOP COMPARISON ***********************/

      *outStream << "\n";
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
