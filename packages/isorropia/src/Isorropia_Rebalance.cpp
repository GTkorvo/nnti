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
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#include <Isorropia_Rebalance.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra_utils.hpp>

#include <Teuchos_RefCountPtr.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#endif

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {

#ifdef HAVE_EPETRA

Epetra_CrsMatrix* create_balanced_copy(const Epetra_CrsMatrix& input_matrix)
{
  //first, create a weights vector which contains the number of nonzeros
  //per row in the input_matrix.
  Epetra_Vector* weights = 0;
  try {
    weights = Isorropia::Epetra_Utils::create_row_weights_nnz(input_matrix);
  }
  catch(std::exception& exc) {
    std::string str1("create_balanced_copy: caught exception: ");
    std::string str2(exc.what());
    throw Isorropia::Exception(str1+str2);
  }

  //now, call the other overloading of 'create_balanced_copy', which
  //accepts a weights vector.
  Epetra_CrsMatrix* balanced_matrix = 0;
  try {
    balanced_matrix =
      Isorropia::create_balanced_copy(input_matrix, *weights);
  }
  catch(std::exception& exc) {
    delete weights;
    throw exc;
  }

  delete weights;

  return balanced_matrix;
}

Epetra_CrsMatrix* create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                                       const Epetra_Vector& row_weights)
{
  //first, figure out what the balanced row-distribution should be, and
  //create a Epetra_Map that describes that row-distribution.
  const Epetra_Comm& comm = input_matrix.Comm();
  //construct a dummy map that will be replaced by the result of the
  //Isorropia function next...
  Epetra_Map bal_rowmap(10, 0, comm);
  try {
    Epetra_Map tmp_map =
      Isorropia::Epetra_Utils::create_rowmap_balanced(input_matrix,
                                                      row_weights);
    bal_rowmap = tmp_map;
  }
  catch(std::exception& exc) {
    std::string str1("create_balanced_copy: caught exception: ");
    std::string str2(exc.what());
    throw Isorropia::Exception(str1+str2);
  }

  //next, create a new Epetra_CrsMatrix (which will be the return-value of
  //this function) with the new row-distribution.
  Epetra_CrsMatrix* balanced_matrix =
    new Epetra_CrsMatrix(Copy, bal_rowmap, 0);

  //finally, import input_matrix into balanced_matrix.
  try {
    Isorropia::Epetra_Utils::import_matrix(input_matrix, *balanced_matrix);
  }
  catch(std::exception& exc) {
    std::string str1("create_balanced_copy: caught exception: ");
    std::string str2(exc.what());
    throw Isorropia::Exception(str1+str2);
  }

  return balanced_matrix;
}

#endif //HAVE_EPETRA

}//namespace Isorropia

