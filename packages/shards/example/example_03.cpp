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
// Questions? Contact Pavel Bochev      (pbboche@sandia.gov)
//                    H. Carter Edwards (hcedwar@sandia.gov)
//                    Denis Ridzal      (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER


/** \file
    \brief  Example of the CellTools class.
    \author Created by P. Bochev, H. Carter Edwards and D. Ridzal
*/
#include <iostream>
#include "Shards_CellTopology.hpp"


using namespace std;
using namespace shards;

/** \brief  Prints the vector with the selected topologies.

    \param  topologies      [in]    - vector containing the selected topologies
    \param  cellType        [in]    - enum for the selected cell type 
    \param  topologyType    [in]    - enum for the selected topology type
  */
void printSelectTopologies(const std::vector<CellTopology>&   topologies,
                     const ECellType                    cellType = ALL_CELLS,
                     const ETopologyType                topologyType = ALL_TOPOLOGIES);




int main(int argc, char *argv[]) {
  std::cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|     Example use of the Shards package:                                      |\n" \
  << "|                                                                             |\n" \
  << "|    1) Query of the available cell topologies                                |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov)                      |\n" \
  << "|                      H. Carter Edwards (hcedwar@sandia.gov)                 |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  shards's website:   http://trilinos.sandia.gov/packages/shards             |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n\n";
  
  
  /*************************************************************************************************
    *                                                                                              *
    *  Cell type enums:     ALL_CELLS,       STANDARD_CELL,  NONSTANDARD_CELL                      *
    *  Topology type enums: ALL_TOPOLOGIES,  BASE_TOPOLOGY,  EXTENDED_TOPOLOGY                     *
    *                                                                                              *
    ***********************************************************************************************/
  std::vector<CellTopology> topologies;
  
  std::cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 1: Queries of predefined 0D, 1D, 2D and 3D cell topologies          |\n"\
    << "===============================================================================\n";
  
  for(int cellDim = 0; cellDim < 4; cellDim++){
    
    std::cout << "************ Selected cell dimension = " << cellDim << " ************\n";
    
    // All cells 
    shards::getTopologies(topologies, cellDim);
    printSelectTopologies(topologies);
    
    
    // All standard cells
    shards::getTopologies(topologies, cellDim, STANDARD_CELL);
    printSelectTopologies(topologies,          STANDARD_CELL);
    
    // All standard cells with base topology
    shards::getTopologies(topologies, cellDim, STANDARD_CELL, BASE_TOPOLOGY);
    printSelectTopologies(topologies,          STANDARD_CELL, BASE_TOPOLOGY);
    
    // All standard cells with extended topology
    shards::getTopologies(topologies, cellDim, STANDARD_CELL, EXTENDED_TOPOLOGY);
    printSelectTopologies(topologies,          STANDARD_CELL, EXTENDED_TOPOLOGY);
    
    
    
    // All non-standard cells
    shards::getTopologies(topologies, cellDim, NONSTANDARD_CELL);
    printSelectTopologies(topologies,          NONSTANDARD_CELL);
    
    // All non-standard 0D cells with base topology
    shards::getTopologies(topologies, cellDim, NONSTANDARD_CELL, BASE_TOPOLOGY);
    printSelectTopologies(topologies,          NONSTANDARD_CELL, BASE_TOPOLOGY);
    
    
    // All non-standard cells with extended topology
    shards::getTopologies(topologies, cellDim,  NONSTANDARD_CELL, EXTENDED_TOPOLOGY);
    printSelectTopologies(topologies,           NONSTANDARD_CELL, EXTENDED_TOPOLOGY);
    
    
  }
  
 
  
  std::cout \
    << "===============================================================================\n"\
    << "| EXAMPLE 2: Query of all predefined cell topologies                          |\n"\
    << "===============================================================================\n";
  
  // This query uses default argument values for all input arguments:
  shards::getTopologies(topologies);
  printSelectTopologies(topologies);
  
  
return 0;
}


/***************************************************************************************************
  *                                                                                                *
  *    Helper function                                                                             *
  *                                                                                                *
  *************************************************************************************************/
void printSelectTopologies(const std::vector<CellTopology>&   topologies,
                     const ECellType                          cellType,
                     const ETopologyType                      topologyType)
{
  
  std::cout << "List of " << shards::ECellTypeToString(cellType) << " ";
  
  // If topologies contains all 33 predefined topologies do not print cell dimension
  if( topologies.size() == 33 ) {
    std::cout << "cells and ";
  }
  else {
    std::cout << topologies[0].getDimension() << "D cells and ";
 
  }
  std::cout << shards::ETopologyTypeToString(topologyType) << " topology types  (total of " 
  << topologies.size() << " cells)\n\n";

  
  for(unsigned i = 0; i < topologies.size(); i++){
    std::cout << topologies[i].getName() << "\n"; 
  }
  std::cout << "===============================================================================\n\n";
}  
























