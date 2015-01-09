/*
//@HEADER
// ***********************************************************************
//
//                     Aristos Optimization Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include <iostream>

/** \class Usr_Par
    \brief Computes and stores several quantities used in the local FE assembly. There is only
    one method \e Print, in addition to the constructor. All computations are performed
    during construction.
*/
            
class Usr_Par {
public:

  Usr_Par();

  Epetra_SerialDenseMatrix Nodes;
  Epetra_SerialDenseVector Weights;

  Epetra_SerialDenseMatrix N;

  Epetra_SerialDenseMatrix Nx1;

  Epetra_SerialDenseMatrix Nx2;

  Epetra_SerialDenseMatrix S1;
  Epetra_SerialDenseMatrix S2;
  Epetra_SerialDenseMatrix S3;

  Epetra_SerialDenseVector Nw;

  Epetra_SerialDenseMatrix NNw;

  Epetra_SerialDenseMatrix * NNNw;

  Epetra_SerialDenseMatrix * NdNdx1Nw;

  Epetra_SerialDenseMatrix * NdNdx2Nw;

  ~Usr_Par() {
    delete [] NNNw;
    delete [] NdNdx1Nw;
    delete [] NdNdx2Nw;
  }
  
  void Print(ostream& os) const;
};
