//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

************************************************************************
*/
//@HEADER

#include <Isorropia_EpetraLibrary.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

#endif

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

namespace Isorropia {

#ifdef HAVE_EPETRA

namespace Epetra {


Library::
Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph, int itype)
  : input_type_(itype),
    numPartSizes(0),
    partGIDs(NULL),
    partSizes(NULL),
    input_graph_(input_graph),
    input_matrix_(0),
    input_coords_(),
    costs_(0),
    weights_(0)
{
  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
}

Library::
Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph, 
	Teuchos::RCP<const Epetra_MultiVector> input_coords, int itype)
  : input_type_(itype),
    numPartSizes(0),
    partGIDs(NULL),
    partSizes(NULL),
    input_graph_(input_graph),
    input_matrix_(0),
    input_coords_(input_coords),
    costs_(0),
    weights_(0)
{
  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
}

Library::
Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
        Teuchos::RCP<CostDescriber> costs, int itype)
  : input_type_(itype),
    numPartSizes(0),
    partGIDs(NULL),
    partSizes(NULL),
    input_graph_(input_graph),
    input_matrix_(0),
    input_coords_(0),
    costs_(costs),
    weights_(0)
{
  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
}

Library::
Library(Teuchos::RCP<const Epetra_CrsGraph> input_graph, Teuchos::RCP<CostDescriber> costs, 
	Teuchos::RCP<const Epetra_MultiVector> input_coords, Teuchos::RCP<const Epetra_MultiVector> weights,
        int itype)
  : input_type_(itype),
    numPartSizes(0),
    partGIDs(NULL),
    partSizes(NULL),
    input_graph_(input_graph),
    input_matrix_(0),
    input_coords_(input_coords),
    costs_(costs),
    weights_(weights)
{
  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
}

Library::
Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, int itype)
  : input_type_(itype),
    numPartSizes(0),
    partGIDs(NULL),
    partSizes(NULL),
    input_graph_(0),
    input_matrix_(input_matrix),
    input_coords_(0),
    costs_(0),
    weights_(0)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
}

Library::
Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, 
	Teuchos::RCP<const Epetra_MultiVector> input_coords,int itype)
  : input_type_(itype),
    numPartSizes(0),
    partGIDs(NULL),
    partSizes(NULL),
    input_graph_(0),
    input_matrix_(input_matrix),
    input_coords_(input_coords),
    costs_(0),
    weights_(0)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
}

Library::
Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	Teuchos::RCP<CostDescriber> costs, int itype)
  : input_type_(itype),
    numPartSizes(0),
    partGIDs(NULL),
    partSizes(NULL),
    input_graph_(0),
    input_matrix_(input_matrix),
    input_coords_(0),
    costs_(costs),
    weights_(0)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
}

////////////////////////////////////////////////////////////////////////////////
Library::
Library(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, Teuchos::RCP<CostDescriber> costs, 
        Teuchos::RCP<const Epetra_MultiVector> input_coords, Teuchos::RCP<const Epetra_MultiVector> weights,
        int itype)
  : input_type_(itype),
    numPartSizes(0),
    partGIDs(NULL),
    partSizes(NULL),
    input_graph_(0),
    input_matrix_(input_matrix),
    input_coords_(input_coords),
    costs_(costs),
    weights_(weights)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
}

Library::
Library(Teuchos::RCP<const Epetra_MultiVector> input_coords, int itype)
  : input_type_(itype),
    numPartSizes(0),
    partGIDs(NULL),
    partSizes(NULL),
    input_graph_(0),
    input_matrix_(0),
    input_coords_(input_coords),
    costs_(0),
    weights_(0)
{
  input_map_ = Teuchos::rcp(&(input_coords->Map()), false);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Library::
Library(Teuchos::RCP<const Epetra_MultiVector> input_coords,
        Teuchos::RCP<const Epetra_MultiVector> weights, int itype)
  : input_type_(itype),
    numPartSizes(0),
    partGIDs(NULL),
    partSizes(NULL),
    input_graph_(0),
    input_matrix_(0),
    input_coords_(input_coords),    
    costs_(0) ,
    weights_(weights)
{
  input_map_ = Teuchos::rcp(&(input_coords->Map()), false);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Library::
Library(Teuchos::RCP<const Epetra_BlockMap> input_map, int itype)
  : input_type_(itype),
    numPartSizes(0),
    partGIDs(NULL),
    partSizes(NULL),
    input_map_(input_map),
    input_graph_(0),
    input_matrix_(0),
    input_coords_(0),
    costs_(0),
    weights_(0)
{

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Library::~Library()
{
  if (partGIDs)
    delete [] partGIDs;
  if (partSizes)
    delete [] partSizes;
}
////////////////////////////////////////////////////////////////////////////////

int Library::precompute()
{
  std::string str1("Isorropia::Epetra::Operator::precompute ");
  std::string str2;

  int inputCount = ((input_graph_.get() == 0) ? 0 : 1);
  inputCount += ((input_matrix_.get() == 0) ? 0 : 1);
  inputCount += ((input_coords_.get() == 0) ? 0 : 1);

  if (inputCount != 1){
    str2 = "ERROR: not holding valid input.";
    throw Isorropia::Exception(str1+str2);
  }

  ///////////////////////////////////
  // If graph info is needed
  ///////////////////////////////////
  if (input_type_ == graph_input_ || input_type_ == hgraph_graph_input_ ||
      input_type_ == graph_geometric_input_ || input_type_ == hgraph_graph_geometric_input_) 
  {
    bool square = false;
    bool symmetric = true;  // no easy way to test for this ?? TODO
    if (input_graph_.get() != 0){
      if (input_graph_->NumGlobalRows() == input_graph_->NumGlobalCols()){
	square = true;
      }
    }
    else if (input_matrix_.get() != 0){
      if (input_matrix_->NumGlobalRows() == input_matrix_->NumGlobalCols()){
	square = true;
      }
    }
    else{
      str2 = "Library requires graph or matrix input";
      throw Isorropia::Exception(str1+str2);
    }
    if (!square){
      str2 = "LB_METHOD=GRAPH, matrix or graph must be square";
      throw Isorropia::Exception(str1+str2);
      return (-1);
    }
    if (!symmetric){
      str2 = "LB_METHOD=GRAPH, matrix or graph must be symmetric";
      throw Isorropia::Exception(str1+str2);
      return (-1);
    }
  }

  ///////////////////////////////////
  // If hypergraph info is needed
  ///////////////////////////////////
  if (input_type_ == hgraph_input_       || input_type_ == hgraph2d_finegrain_input_ ||
      input_type_ == hgraph_graph_input_ || input_type_ == hgraph_geometric_input_ ||
      input_type_ == hgraph_graph_geometric_input_)
  {
    if ((input_graph_.get() == 0) && (input_matrix_.get() == 0))
    {
      str2 = "Library requires graph or matrix input";
      throw Isorropia::Exception(str1+str2);
    }
  }

  ///////////////////////////////////
  // If geometric info is needed
  ///////////////////////////////////
  if (input_type_ == geometric_input_ || input_type_ == hgraph_geometric_input_ ||
      input_type_ == graph_geometric_input_ || input_type_ == hgraph_graph_geometric_input_)
  {
    if ((input_coords_.get() == 0) ||
        (input_coords_->NumVectors() < 1) || (input_coords_->NumVectors() > 3)){
      str2 = "Operation requires 1, 2 or 3 dimensional coordinate input";
      throw Isorropia::Exception(str1+str2);
    }

    if (weights_.get() != 0){
      if (weights_->MyLength() != input_coords_->MyLength()){
        str2 = "Number of weights does not equal number of coordinates";
        throw Isorropia::Exception(str1+str2);
      }
    }
  }
  return (0);
}


} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

