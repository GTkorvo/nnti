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

#ifndef _Isorropia_Rebalance_hpp_
#define _Isorropia_Rebalance_hpp_

#include <Isorropia_configdefs.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_LinearProblem;
#endif

/** Isorropia is the namespace that contains isorropia's declarations
  for classes and functions.
*/
namespace Isorropia {

#ifdef HAVE_EPETRA
/** Create a balanced copy of an input Epetra_CrsMatrix.

  This function represents the most basic default case, not accepting
  parameters or weights, and uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
  create_balanced_copy(const Epetra_CrsMatrix& input_matrix);

/** Create a balanced copy of an input Epetra_RowMatrix.

  This function represents the most basic default case, not accepting
  parameters or weights, and uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
  create_balanced_copy(const Epetra_RowMatrix& input_matrix);

/** Create a balanced copy of an input Epetra_CrsMatrix, accounting for
   user-supplied weights assigned to each row.

  This function uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to assign 'chunks' of rows
  such that the sum of associated weights is equal on each partition.

  The row_weights vector is required to have the same size and distribution
  as the row-map of input_matrix.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
  create_balanced_copy(const Epetra_CrsMatrix& input_matrix,
                       const Epetra_Vector& row_weights);

/** Create a balanced copy of an input Epetra_RowMatrix, accounting for
   user-supplied weights assigned to each row.

  This function uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to assign 'chunks' of rows
  such that the sum of associated weights is equal on each partition.

  The row_weights vector is required to have the same size and distribution
  as the row-map of input_matrix.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
  create_balanced_copy(const Epetra_RowMatrix& input_matrix,
                       const Epetra_Vector& row_weights);

/** Create a balanced copy of an input Epetra_CrsGraph.

  This function represents a basic case, not accepting weights.
  By default, it uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  If the ParameterList object contains a string parameter with
  name == "Balancing package" and value == "Zoltan", then the Zoltan
  library is used to perform the balancing.

  Furthmore, if the ParameterList object contains a string parameter with
  name == "Partitioning algorithm" then possible values are "graph"
  or "hypergraph" (not yet supported). If "Zoltan" has not been specified
  as the "Balancing package", then "Partitioning algorithm" is ignored.

  The rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros equal in each partition. I.e., it is equivalent to specifying
  a weighted rebalance where the weights are the number of nonzeros in
  each row.
*/
Teuchos::RefCountPtr<Epetra_CrsGraph>
create_balanced_copy(const Epetra_CrsGraph& input_graph,
		     Teuchos::ParameterList& paramlist);

/** Create a balanced copy of an input Epetra_CrsGraph, accounting for
   user-supplied weights assigned to each row.

  This function uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to assign 'chunks' of rows
  such that the sum of associated weights is equal on each partition.

  The row_weights vector is required to have the same size and distribution
  as the row-map of input_graph.
*/

Teuchos::RefCountPtr<Epetra_CrsGraph>
  create_balanced_copy(const Epetra_CrsGraph& input_graph,
                       const Epetra_Vector& row_weights);

/** Create a balanced copy of an input Epetra_LinearProblem.

  This function represents the most basic default case, not accepting
  parameters or weights, and uses a self-contained internal implementation
  rather than interfacing to a third-party library such as Zoltan.

  The rebalancing is 1-D, row-wise, and attempts to make the number
  of nonzeros in the matrix equal in each partition. I.e., it is equivalent
  to specifying a weighted matrix rebalance where the weights are the number
  of nonzeros in each row. The vectors are then redistributed to match the
  row-distribution of the matrix.

  Important note: It is the caller's responsibility to destroy the
  matrix and vector attributes of the returned Epetra_LinearProblem.
*/
Teuchos::RefCountPtr<Epetra_LinearProblem>
  create_balanced_copy(const Epetra_LinearProblem& input_problem);

/** Return a new Epetra_CrsMatrix object constructed with target_rowmap,
  and with the contents of input_matrix imported into it.

  The caller is responsible for deleting the returned object.

  @param input_matrix Source/input object.

  @param target_rowmap Target rowmap, required to be compatible with
     input_matrix.RowMap() in terms of number-of-elements, etc.

  @param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
  redistribute_rows(const Epetra_CrsMatrix& input_matrix,
                    const Epetra_Map& target_rowmap,
                    Epetra_Import* importer=0);

/** Return a new Epetra_CrsMatrix object constructed with target_rowmap,
  and with the contents of input_matrix imported into it.

  The caller is responsible for deleting the returned object.

  @param input_matrix Source/input object.

  @param target_rowmap Target rowmap, required to be compatible with
     input_matrix.RowMatrixRowMap() in terms of number-of-elements, etc.

  @param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RefCountPtr<Epetra_CrsMatrix>
  redistribute_rows(const Epetra_RowMatrix& input_matrix,
                    const Epetra_Map& target_rowmap,
                    Epetra_Import* importer=0);

/** Return a new Epetra_CrsGraph object constructed with target_rowmap,
  and with the contents of input_graph imported into it.

  @param input_graph Source/input object.

  @param target_rowmap Target rowmap, required to be compatible with
     input_graph.RowMap() in terms of number-of-elements, etc.

  @param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RefCountPtr<Epetra_CrsGraph>
  redistribute_rows(const Epetra_CrsGraph& input_graph,
                    const Epetra_Map& target_rowmap,
                    Epetra_Import* importer=0);

/** Return a new Epetra_MultiVector object constructed with target_map,
  and with the contents of 'input' imported into it.

  @param input Source/input object.

  @param target_map Target map, required to be compatible with
     input.Map() in terms of number-of-elements, etc.

  @param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RefCountPtr<Epetra_MultiVector>
  redistribute(const Epetra_MultiVector& input,
               const Epetra_BlockMap& target_map,
               Epetra_Import* importer=0);

/** Return a new Epetra_Vector object constructed with target_map,
  and with the contents of 'input' imported into it.

  @param input Source/input object.

  @param target_map Target map, required to be compatible with
     input.RowMap() in terms of number-of-elements, etc.

  @param importer Optional argument. If importer is supplied, it will be
     used to perform the import operation. Otherwise, a temporary importer
     will be created and used.
*/
Teuchos::RefCountPtr<Epetra_Vector>
  redistribute(const Epetra_Vector& input,
               const Epetra_Map& target_map,
               Epetra_Import* importer=0);

#endif //HAVE_EPETRA

}//namespace Isorropia

#endif

