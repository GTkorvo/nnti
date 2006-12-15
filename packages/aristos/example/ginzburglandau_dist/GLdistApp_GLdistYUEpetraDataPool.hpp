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
  Teuchos::RefCountPtr<Epetra_SerialDenseMatrix> ipcoords_;
  /** \brief Global nodes (interior, nonoverlapping) in this subdomain.*/
  Teuchos::RefCountPtr<Epetra_IntSerialDenseVector> ipindx_;
  /** \brief Coordinates of all nodes in this subdomain.*/
  Teuchos::RefCountPtr<Epetra_SerialDenseMatrix> pcoords_;
  /** \brief Global nodes (interior + shared, overlapping) in this subdomain.*/
  Teuchos::RefCountPtr<Epetra_IntSerialDenseVector> pindx_;
  /** \brief Elements (this includes all overlapping nodes).*/
  Teuchos::RefCountPtr<Epetra_IntSerialDenseMatrix> t_;
  /** \brief Edges.*/
  Teuchos::RefCountPtr<Epetra_IntSerialDenseMatrix> e_;

  /** \brief Volume stiffness matrix.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> A_;
  /** \brief Control/state mass matrix.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> B_;
  /** \brief Volume mass matrix.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> H_;
  /** \brief Edge mass matrix.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> R_;

  /** \brief Augmented system matrix: \n
   [ I  Jac* ] \n
   [Jac  0   ]
  */
  Teuchos::RefCountPtr<Epetra_CrsMatrix> Augmat_;

  /** \brief DD preconditioning operator.*/
  //Ifpack_Preconditioner* Prec_;
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Prec_;

  /** \brief Jacobian of the nonlinear term.*/
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> Npy_;

  /** \brief Right-hand side of the PDE.*/
  Teuchos::RefCountPtr<Epetra_FEVector> b_;
  /** \brief The desired state.*/
  Teuchos::RefCountPtr<Epetra_FEVector> q_;

  Teuchos::RefCountPtr<Epetra_FEVector> Ny_;

  /** \brief Regularization parameter.*/
  double beta_;

  /** \brief Base geometry file name.*/
  char geomfile_[120];

public:

  GLdistYUEpetraDataPool( Epetra_Comm * commptr, double beta, char * myfile );

  /** \brief Calls functions to compute nonlinear quantities and the augmented system matrix.
             These computations are performed after every update of the SQP iterate.
  */
  void computeAll( const Aristos::Vector &x );

  /** \brief Solves augmented system.*/
  int  solveAugsys( const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsy,
                    const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsu,
                    const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsp,
                    const Teuchos::RefCountPtr<Epetra_MultiVector> & y,
                    const Teuchos::RefCountPtr<Epetra_MultiVector> & u,
                    const Teuchos::RefCountPtr<Epetra_MultiVector> & p,
                    double * tol );

  /** \brief Solves augmented system, with dynamic tolerances.*/
  int  solveAugsysDyn( const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsy,
                       const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsu,
                       const Teuchos::RefCountPtr<const Epetra_MultiVector> & rhsp,
                       const Teuchos::RefCountPtr<Epetra_MultiVector> & y,
                       const Teuchos::RefCountPtr<Epetra_MultiVector> & u,
                       const Teuchos::RefCountPtr<Epetra_MultiVector> & p,
                       double * tol );

  Epetra_Comm* getCommPtr();

  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getA();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getB();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getH();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getR();
  Teuchos::RefCountPtr<Epetra_CrsMatrix> getAugmat();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> getNpy();

  Teuchos::RefCountPtr<Epetra_FEVector> getb();
  Teuchos::RefCountPtr<Epetra_FEVector> getq();
  Teuchos::RefCountPtr<Epetra_FEVector> getNy();

  /** \brief Calls the function that computes the nonlinear term.*/
  void computeNy(const Teuchos::RefCountPtr<const Epetra_MultiVector> & y);

  /** \brief Calls the function that computes the Jacobian of the nonlinear term.*/
  void computeNpy(const Teuchos::RefCountPtr<const Epetra_MultiVector> & y);

  /** \brief Assembles the augmented system (KKT-type) matrix.*/
  void computeAugmat();

  /** \brief Computes DD Ifpack preconditioner.*/
  int computePrec();
  
  Teuchos::RefCountPtr<const Epetra_SerialDenseMatrix> getipcoords();
  Teuchos::RefCountPtr<const Epetra_IntSerialDenseVector> getipindx();
  Teuchos::RefCountPtr<const Epetra_SerialDenseMatrix> getpcoords();
  Teuchos::RefCountPtr<const Epetra_IntSerialDenseVector> getpindx();
  Teuchos::RefCountPtr<const Epetra_IntSerialDenseMatrix> gett();
  Teuchos::RefCountPtr<const Epetra_IntSerialDenseMatrix> gete();

  double getbeta();
  
  /** \brief Outputs the solution vector to Matlab files.*/
  void PrintSolutionMatlab( const Teuchos::RefCountPtr<const Epetra_Vector> & x );

  /** \brief Outputs the solution vector to VTK files.*/
  void PrintSolutionVTK( const Teuchos::RefCountPtr<const Epetra_Vector> & y );
};

} // namespace GLdistApp

#endif
