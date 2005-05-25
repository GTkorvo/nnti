// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Amesos_config.h"
#include "Amesos.h"
#include "Amesos_Klu.h"
#ifdef HAVE_AMESOS_LAPACK
#include "Amesos_Lapack.h"
#endif
#ifdef HAVE_AMESOS_MUMPS
#include "Amesos_Mumps.h"
#endif
#ifdef HAVE_AMESOS_SCALAPACK
#include "Amesos_Scalapack.h"
#endif
#ifdef HAVE_AMESOS_UMFPACK
#include "Amesos_Umfpack.h"
#endif
#ifdef HAVE_AMESOS_SUPERLUDIST
#include "Amesos_Superludist.h"
#endif
#ifdef HAVE_AMESOS_SUPERLU
#include "Amesos_Superlu.h"
#endif
#ifdef HAVE_AMESOS_DSCPACK
#include "Amesos_Dscpack.h"
#endif
#ifdef HAVE_AMESOS_PARDISO
#include "Amesos_Pardiso.h"
#endif
#ifdef HAVE_AMESOS_TAUCS
#include "Amesos_Taucs.h"
#endif
#include "Epetra_Object.h"

static bool verbose = false; 

Amesos_BaseSolver* Amesos::Create(const char* ClassType, 
				  const Epetra_LinearProblem& LinearProblem ) 
{ 
  string CT = ClassType; 
  return(Create(CT,LinearProblem));
}

Amesos_BaseSolver* Amesos::Create(const string CT,
				  const Epetra_LinearProblem& LinearProblem )
{

  if ((CT == "Amesos_Lapack") || (CT == "Lapack")) { 
#ifdef HAVE_AMESOS_LAPACK
    return new Amesos_Lapack(LinearProblem); 
#else
    if (verbose) cerr << "Amesos_Lapack is not implemented" << endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Klu") || (CT == "Klu")) { 
#ifdef HAVE_AMESOS_KLU
    return new Amesos_Klu(LinearProblem); 
#else
    if (verbose) cerr << "Amesos_Klu is not implemented" << endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Umfpack") || (CT == "UmfpacK")) { 
#ifdef HAVE_AMESOS_UMFPACK
    return new Amesos_Umfpack(LinearProblem); 
#else
    if (verbose) cerr << "Amesos_Umfpack is not implemented" << endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Superlu") || (CT == "Superlu")) { 
#ifdef HAVE_AMESOS_SUPERLU
    return new Amesos_Superlu(LinearProblem); 
#else
    if (verbose) cerr << "Amesos_Superlu is not implemented" << endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Superludist") || (CT == "Superludist")) { 
#ifdef HAVE_AMESOS_SUPERLUDIST
    return new Amesos_Superludist(LinearProblem); 
#else
    if (verbose) cerr << "Amesos_Superludist is not implemented" << endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Mumps") || (CT == "Mumps")) { 
#ifdef HAVE_AMESOS_MUMPS
    return new Amesos_Mumps(LinearProblem); 
#else
    if (verbose) cerr << "Amesos_Mumps is not implemented" << endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Scalapack") || (CT == "Scalapack")) { 
#ifdef HAVE_AMESOS_SCALAPACK
    return new Amesos_Scalapack(LinearProblem); 
#else
    if (verbose) cerr << "Amesos_Scalapack is not implemented" << endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Dscpack") || (CT == "Dscpack")) { 
#ifdef HAVE_AMESOS_DSCPACK
    return new Amesos_Dscpack(LinearProblem); 
#else
    if (verbose) cerr << "Amesos_Dscpack is not implemented" << endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Pardiso") || (CT == "Pardiso")) { 
#ifdef HAVE_AMESOS_PARDISO
    return new Amesos_Pardiso(LinearProblem); 
#else
    if (verbose) cerr << "Amesos_Pardiso is not implemented" << endl ; 
    return(0); 
#endif
  } 
  
  if ((CT == "Amesos_Taucs") || (CT == "Taucs")) { 
#ifdef HAVE_AMESOS_TAUCS
    return new Amesos_Taucs(LinearProblem); 
#else
    if (verbose) cerr << "Amesos_Taucs is not implemented" << endl ; 
    return(0); 
#endif
  } 
  
  if (verbose) cerr << "Unknown class type:" << CT << endl ; 
  return(0); 
}

// ====================================================================
bool Amesos::Query(const char* ClassType)
{
  string CT = ClassType;
  return(Query(CT));
}

// ====================================================================
bool Amesos::Query(const string CT) 
{ 

  if ((CT == "Amesos_Lapack") || (CT == "Lapack")) { 
#ifdef HAVE_AMESOS_LAPACK
    return true;
#else
    return false;
#endif
  } 
  
  if ((CT == "Amesos_Klu") || (CT == "Klu")) { 
#ifdef HAVE_AMESOS_KLU
    return true; 
#else
    return false;
#endif
  } 
  
  if ((CT == "Amesos_Umfpack") || (CT == "Umfpack")) { 
#ifdef HAVE_AMESOS_UMFPACK
    return true;
#else
    return false;
#endif
  } 
  
  if ((CT == "Amesos_Superlu") || ( CT == "Superlu")) { 
#ifdef HAVE_AMESOS_SUPERLU
    return true; 
#else
    return false;
#endif
  }

  if ((CT == "Amesos_Superludist") || (CT == "Superludist")) { 
#ifdef HAVE_AMESOS_SUPERLUDIST
    return true;
#else
    return false;
#endif
  } 

  if ((CT == "Amesos_Mumps") || (CT == "Mumps")) { 
#ifdef HAVE_AMESOS_MUMPS
    return true;
#else
    return false;
#endif
  } 

  if ((CT == "Amesos_Scalapack") || (CT == "Scalapack")) { 
#ifdef HAVE_AMESOS_SCALAPACK
    return true;
#else
    return false;
#endif
  } 

  if ((CT == "Amesos_Dscpack") || (CT == "Dscpack")) { 
#ifdef HAVE_AMESOS_DSCPACK
    return true;
#else
    return false;
#endif
  } 
  
  if ((CT == "Amesos_Pardiso") || (CT == "Pardiso")) { 
#ifdef HAVE_AMESOS_PARDISO
    return true;
#else
    return false;
#endif
  } 
  
  if ((CT == "Amesos_Taucs") || (CT == "Taucs")) { 
#ifdef HAVE_AMESOS_TAUCS
    return true;
#else
    return false;
#endif
  } 
  
  return(false);

}
