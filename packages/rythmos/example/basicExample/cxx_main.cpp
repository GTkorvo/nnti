//@HEADER
// ************************************************************************
// 
//                          Rythmos Package 
//                 Copyright (2005) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"

#include "ExampleApplication.hpp"

// Includes for Rythmos:
#include "Rythmos_ConfigDefs.h"
//#include "ExampleApplicationRythmosInterface.hpp"
#include "Rythmos_ForwardEulerStepper.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"

// Includes for Thyra:
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"

#include <string>

// Includes for Teuchos:
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

int main(int argc, char *argv[])
{
  bool verbose = true; // verbosity level.
  bool success = true; // determine if the run was successfull

  try { // catch exceptions

#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
    MPI_Comm mpiComm = MPI_COMM_WORLD;
    int procRank = 0;
    int numProc;
    MPI_Comm_size( mpiComm, &numProc );
    MPI_Comm_rank( mpiComm, &procRank );
#endif // HAVE_MPI

    double lambda_min = -0.9;   // min ODE coefficient
    double lambda_max = -0.01;  // max ODE coefficient
    std::string lambda_fit = "random"; // Lambda model
    int numElements = 1; // number of elements in vector
    double x0 = 10.0; // ODE initial condition
    int N = 10;  // number of steps to take
    std::string method = "FE";  // other choice is method="ERK4"
    bool version = false;  // display version information 

    // Parse the command-line options:
    Teuchos::CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption( "x0", &x0, "Constant ODE initial condition." );
    clp.setOption( "lambda_min", &lambda_min, "Lower bound for ODE coefficient");
    clp.setOption( "lambda_max", &lambda_max, "Upper bound for ODE coefficient");
    clp.setOption( "lambda_fit", &lambda_fit, "Lambda model:  random, linear");
    clp.setOption( "numelements", &numElements, "Problem size");
    clp.setOption( "method", &method, "Integration method:  FE, ERK4." );
    clp.setOption( "numsteps", &N, "Number of integration steps to take" );
    clp.setOption( "verbose", "quiet", &verbose, "Set if output is printed or not" );
    clp.setOption( "version", "run", &version, "Version of this code" );

    Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    if (version) // Display version information and exit.
    {
      std::cout << Rythmos::Rythmos_Version() << std::endl; 
      std::cout << "basicExample Version 0.1 - 06/23/05" << std::endl;
      return(0);
    }

    if (lambda_min > lambda_max)
    {
      std::cerr << "lamba_min must be less than lambda_max" << std::endl;
      return(1);
    }
    
    // Set up the parameter list for the application:
    Teuchos::ParameterList params;
    params.set( "Lambda_min", lambda_min );
    params.set( "Lambda_max", lambda_max );
    params.set( "Lambda_fit", lambda_fit );
    params.set( "NumElements", numElements );
    params.set( "x0", x0 );
#ifdef HAVE_MPI
    params.set( "MPIComm", mpiComm );
#endif // HAVE_MPI

    // create interface to problem
    Teuchos::RefCountPtr<ExampleApplication>
      epetraModel = Teuchos::rcp(new ExampleApplication(params));
    Teuchos::RefCountPtr<Thyra::ModelEvaluator<double> >
      model = Teuchos::rcp(new Thyra::EpetraModelEvaluator(epetraModel));
    
    // Create Stepper object depending on command-line input
    Teuchos::RefCountPtr<Rythmos::Stepper<double> > stepper_ptr;
    if (method == "ERK4") // Explicit Runge-Kutta 4 stage
    {
      stepper_ptr = Teuchos::rcp(new Rythmos::ExplicitRKStepper<double>(model));
      method = "Explicit Runge-Kutta of order 4";
    } else // Forward Euler
    {
      stepper_ptr = Teuchos::rcp(new Rythmos::ForwardEulerStepper<double>(model));
      method = "Forward Euler";
    }
    Rythmos::Stepper<double> &stepper = *stepper_ptr;


    double t0 = 0.0;
    double t1 = 1.0;
    double dt = (t1-t0)/N;

    // Integrate forward with fixed step sizes:
    for (int i=1 ; i<=N ; ++i)
    {
      double dt_taken = stepper.TakeStep(dt);
      if (dt_taken != dt)
      {
        cerr << "Error, stepper took step of dt = " << dt_taken << " when asked to take step of dt = " << dt << std::endl;
        break;
      }
    }
    // Get solution out of stepper:
    Teuchos::RefCountPtr<const Thyra::VectorBase<double> > x_computed_thyra_ptr = stepper.get_solution();
    // Convert Thyra::VectorBase to Epetra_Vector
    Teuchos::RefCountPtr<const Epetra_Vector>
      x_computed_ptr = Thyra::get_Epetra_Vector(*(epetraModel->get_x_map()),x_computed_thyra_ptr);
    const Epetra_Vector &x_computed = *x_computed_ptr;

    // compute exact answer
    Teuchos::RefCountPtr<const Epetra_Vector> lambda_ptr = epetraModel->get_coeff();
    const Epetra_Vector &lambda = *lambda_ptr;
    Epetra_Vector x_star(lambda.Map());
    for (int i=0 ; i < x_star.MyLength() ; ++i)
    {
      x_star[i] = x0*exp(lambda[i]*t1);
    }

    // 06/03/05 tscoffe to get an Epetra_Map associated with an Epetra_Vector:
    // x.Map()
    // to get an Epetra_Comm associated with an Epetra_Vector:
    // x.Comm()
    
  //  Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm = (*epetraModel).get_epetra_comm();
    //int MyPID = epetraModel->get_Epetra_Map()->Comm()->MyPID();
    int MyPID = x_computed.Comm().MyPID();
  //  int MyPID = epetra_comm->MyPID();
    if (MyPID == 0)
    {
      std::cout << "Integrating \\dot{x}=\\lambda x from t = " << t0 
                << " to t = " << t1 << std::endl;
      std::cout << "using " << method << std::endl;
      std::cout << "with initial x_0 = " << x0
                << ", \\Delta t = " << dt  << "." << std::endl;
    }
    int MyLength = x_computed.MyLength();
    for (int i=0 ; i<MyLength ; ++i)
    {
      std::cout << "lambda[" << MyPID*MyLength+i << "] = " << lambda[i] << std::endl;
    }
    for (int i=0 ; i<MyLength ; ++i)
    {
      std::cout << "Computed: x[" << MyPID*MyLength+i << "] = " << x_computed[i] << "\t" 
                <<    "Exact: x[" << MyPID*MyLength+i << "] = " << x_star[i] << std::endl;
    }
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif // HAVE_MPI

   } // end try
   catch( const std::exception &excpt ) {
    if(verbose)
      std::cerr << "*** Caught a standard exception : " << excpt.what() << std::endl;
    success = false;
  }
  catch( ... ) {
    if(verbose)
      std::cerr << "*** Caught an unknown exception!\n";
    success = false;
  }

  return success ? 0 : 1;
} // end main() [Doxygen looks for this!]

