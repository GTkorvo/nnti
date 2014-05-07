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

#ifndef GLDISTAPP_YUEPETRADATAPOOL_H
#define GLDISTAPP_YUEPETRADATAPOOL_H

//#include "Epetra_config.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "AztecOO_Operator.h"
#include "Epetra_LAPACK.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"

#include "Amesos.h"
#include "Ifpack.h"
#include "Amesos_Klu.h"
#include "Amesos_Umfpack.h"
//#include "Amesos_Superludist.h"
#include "Teuchos_ParameterList.hpp"

#include "usr_par.h"


#include "Aristos_DataPool.hpp"
#include "Aristos_YUEpetraVector.hpp"
#include "includes.h"

namespace GLdistApp {

/**
    \brief Implements the Aristos::DataPool interface module for the parallel
    Ginzburg-Landau (GLdist) application.
*/
class GLdistYUEpetraDataPool : public Aristos::DataPool {
private:

  /** \brief Communicator member.*/
  Epetra_Comm * commptr_;
          
  /** \brief Coordinates of nodes that are unique to this subdomain.*/
  Teuchos::RCP<Epetra_SerialDenseMatrix> ipcoords_;
  /** \brief Global nodes (interior, nonoverlapping) in this subdomain.*/
  Teuchos::RCP<Epetra_IntSerialDenseVector> ipindx_;
  /** \brief Coordinates of all nodes in this subdomain.*/
  Teuchos::RCP<Epetra_SerialDenseMatrix> pcoords_;
  /** \brief Global nodes (interior + shared, overlapping) in this subdomain.*/
  Teuchos::RCP<Epetra_IntSerialDenseVector> pindx_;
  /** \brief Elements (this includes all overlapping nodes).*/
  Teuchos::RCP<Epetra_IntSerialDenseMatrix> t_;
  /** \brief Edges.*/
  Teuchos::RCP<Epetra_IntSerialDenseMatrix> e_;

  /** \brief Volume stiffness matrix.*/
  Teuchos::RCP<Epetra_FECrsMatrix> A_;
  /** \brief Control/state mass matrix.*/
  Teuchos::RCP<Epetra_FECrsMatrix> B_;
  /** \brief Volume mass matrix.*/
  Teuchos::RCP<Epetra_FECrsMatrix> H_;
  /** \brief Edge mass matrix.*/
  Teuchos::RCP<Epetra_FECrsMatrix> R_;

  /** \brief Augmented system matrix: \n
   [ I  Jac* ] \n
   [Jac  0   ]
  */
  Teuchos::RCP<Epetra_CrsMatrix> Augmat_;
  Teuchos::RCP<const Thyra::LinearOpBase<double> > thyra_Augmat_;

  /** \brief DD preconditioning operator.*/
  //Ifpack_Preconditioner* Prec_;
  Teuchos::RCP<Ifpack_Preconditioner> Prec_;

  /** \brief Jacobian of the nonlinear term.*/
  Teuchos::RCP<Epetra_FECrsMatrix> Npy_;

  /** \brief LOWSF object for the Augmat_. */
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > Augmat_lowsFactory_;

  /** \brief LOWS object for the Augmat_. */
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > Augmat_lows_;

  /** \brief Right-hand side of the PDE.*/
  Teuchos::RCP<Epetra_FEVector> b_;
  /** \brief The desired state.*/
  Teuchos::RCP<Epetra_FEVector> q_;

  Teuchos::RCP<Epetra_FEVector> Ny_;

  /** \brief Regularization parameter.*/
  double beta_;

  /** \brief Base geometry file name.*/
  char geomfile_[120];

public:

  // ToDo: Make these private!
  static bool useStratimikos;
  static std::string stratimikosXmlFile;

  GLdistYUEpetraDataPool( Epetra_Comm * commptr, double beta, const std::string &myfile );

  /** \brief Calls functions to compute nonlinear quantities and the augmented system matrix.
             These computations are performed after every update of the SQP iterate.
  */
  void computeAll( const Aristos::Vector &x );

  /** \brief Solves augmented system.*/
  int  solveAugsys( const Teuchos::RCP<const Epetra_MultiVector> & rhsy,
                    const Teuchos::RCP<const Epetra_MultiVector> & rhsu,
                    const Teuchos::RCP<const Epetra_MultiVector> & rhsp,
                    const Teuchos::RCP<Epetra_MultiVector> & y,
                    const Teuchos::RCP<Epetra_MultiVector> & u,
                    const Teuchos::RCP<Epetra_MultiVector> & p,
                    double * tol );

  /** \brief Solves augmented system, with dynamic tolerances.*/
  int  solveAugsysDyn( const Teuchos::RCP<const Epetra_MultiVector> & rhsy,
                       const Teuchos::RCP<const Epetra_MultiVector> & rhsu,
                       const Teuchos::RCP<const Epetra_MultiVector> & rhsp,
                       const Teuchos::RCP<Epetra_MultiVector> & y,
                       const Teuchos::RCP<Epetra_MultiVector> & u,
                       const Teuchos::RCP<Epetra_MultiVector> & p,
                       double * tol );

  Epetra_Comm* getCommPtr();

  Teuchos::RCP<Epetra_FECrsMatrix> getA();
  Teuchos::RCP<Epetra_FECrsMatrix> getB();
  Teuchos::RCP<Epetra_FECrsMatrix> getH();
  Teuchos::RCP<Epetra_FECrsMatrix> getR();
  Teuchos::RCP<Epetra_CrsMatrix> getAugmat();
  Teuchos::RCP<Epetra_FECrsMatrix> getNpy();

  Teuchos::RCP<Epetra_FEVector> getb();
  Teuchos::RCP<Epetra_FEVector> getq();
  Teuchos::RCP<Epetra_FEVector> getNy();

  /** \brief Calls the function that computes the nonlinear term.*/
  void computeNy(const Teuchos::RCP<const Epetra_MultiVector> & y);

  /** \brief Calls the function that computes the Jacobian of the nonlinear term.*/
  void computeNpy(const Teuchos::RCP<const Epetra_MultiVector> & y);

  /** \brief Assembles the augmented system (KKT-type) matrix.*/
  void computeAugmat();

  /** \brief Computes DD Ifpack preconditioner.*/
  int computePrec();
  
  Teuchos::RCP<const Epetra_SerialDenseMatrix> getipcoords();
  Teuchos::RCP<const Epetra_IntSerialDenseVector> getipindx();
  Teuchos::RCP<const Epetra_SerialDenseMatrix> getpcoords();
  Teuchos::RCP<const Epetra_IntSerialDenseVector> getpindx();
  Teuchos::RCP<const Epetra_IntSerialDenseMatrix> gett();
  Teuchos::RCP<const Epetra_IntSerialDenseMatrix> gete();

  double getbeta();
  
  /** \brief Outputs the solution vector to Matlab files.*/
  void PrintSolutionMatlab( const Teuchos::RCP<const Epetra_Vector> & x );

  /** \brief Outputs the solution vector to VTK files.*/
  void PrintSolutionVTK( const Teuchos::RCP<const Epetra_Vector> & y );
};

} // namespace GLdistApp

#endif
