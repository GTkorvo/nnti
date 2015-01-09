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

class Epetra_Comm;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_FECrsMatrix;
class Epetra_FEVector;
class Epetra_SerialDenseMatrix;
class Epetra_SerialDenseVector;
class Epetra_IntSerialDenseVector;
class Epetra_IntSerialDenseMatrix;
namespace Teuchos {
 template <class T> class RefCountPtr;
}
class Usr_Par;


bool CrsMatrix2MATLAB(const Epetra_CrsMatrix &, ostream &);

bool Vector2MATLAB( const Epetra_Vector &, ostream &);

bool FEVector2MATLAB( const Epetra_FEVector &, ostream &);

int quadrature(const int, const int, Epetra_SerialDenseMatrix &,
               Epetra_SerialDenseVector &);

int meshreader(const Epetra_Comm &,
               Epetra_IntSerialDenseVector &,
               Epetra_SerialDenseMatrix &,
               Epetra_IntSerialDenseVector &,
               Epetra_SerialDenseMatrix &,
               Epetra_IntSerialDenseMatrix &,
               Epetra_IntSerialDenseMatrix &,
               const char *);

int lassembly(const Epetra_SerialDenseMatrix &,
              const Epetra_SerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_SerialDenseVector &,
              const Epetra_SerialDenseVector &,
              const Usr_Par &,
              Epetra_SerialDenseMatrix &,
              Epetra_SerialDenseVector &);

int assemblyFECrs(const Epetra_Comm &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
                  Teuchos::RefCountPtr<Epetra_FEVector> &);

int assemblyFECrs(const Epetra_Comm &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
                  Teuchos::RefCountPtr<Epetra_FEVector> &,
                  bool);

int assemble(const Epetra_Comm &,
             const Epetra_IntSerialDenseVector &,
             const Epetra_SerialDenseMatrix &,
             const Epetra_IntSerialDenseVector &,
             const Epetra_SerialDenseMatrix &,
             const Epetra_IntSerialDenseMatrix &,
             const Epetra_IntSerialDenseMatrix &,
             Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
             Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
             Teuchos::RefCountPtr<Epetra_FEVector> &);

int assemble_bdry(const Epetra_Comm &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  Teuchos::RefCountPtr<Epetra_FECrsMatrix> &,
                  Teuchos::RefCountPtr<Epetra_FECrsMatrix> &);

int nonlinvec(const Epetra_Comm &,
              const Epetra_IntSerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_IntSerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_IntSerialDenseMatrix &,
              const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
              Teuchos::RefCountPtr<Epetra_FEVector> &);

int nonlinjac(const Epetra_Comm &,
              const Epetra_IntSerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_IntSerialDenseVector &,
              const Epetra_SerialDenseMatrix &,
              const Epetra_IntSerialDenseMatrix &,
              const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
              Teuchos::RefCountPtr<Epetra_FECrsMatrix> &);

int nonlinhessvec(const Epetra_Comm &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseVector &,
                  const Epetra_SerialDenseMatrix &,
                  const Epetra_IntSerialDenseMatrix &,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
                  const Teuchos::RefCountPtr<const Epetra_MultiVector> &,
                  Teuchos::RefCountPtr<Epetra_FEVector> &);


ostream& operator <<(ostream &, const Usr_Par &);
