/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef PLAYA_AZTECSOLVER_HPP
#define PLAYA_AZTECSOLVER_HPP

#include "PlayaDefs.hpp"
#include "PlayaLinearSolverBaseDecl.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaPrintable.hpp"
#include "Teuchos_Describable.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include <map>

#include "AztecOO.h"

#define AZ_recursive_iterate 10001

#define HAVE_ML


namespace Playa
{
  using namespace Teuchos;

  /**
   *
   */
  class AztecSolver : public LinearSolverBase<double>,
                      public Playa::Handleable<LinearSolverBase<double> >,
                      public Printable,
                      public Describable
  {
  public:
    /** */
    AztecSolver(const Teuchos::map<int, int>& aztecOptions,
                const Teuchos::map<int, double>& aztecParameters);

    /** */
    AztecSolver(const Teuchos::ParameterList& params);

    /** */
    virtual ~AztecSolver(){;}

    /** Change the convergence tolerance. */
    virtual void updateTolerance(const double& tol);


    /** Set the preconditioning operator */
    void setUserPrec(const LinearOperator<double>& P,
		     const LinearSolver<double>& pSolver);

    /** \name Printable interface */
    //@{
    /** Write to a stream  */
    void print(std::ostream& os) const 
    {
      os << description() << std::endl;
    }
    //@}
    
    /** \name Describable interface */
    //@{
    /** Write a brief description */
    std::string description() const {return "AztecSolver";}
    //@}

    

    /** */
    virtual SolverState<double> solve(const LinearOperator<double>& op,
                                      const Vector<double>& rhs,
                                      Vector<double>& soln) const ;

    /** \name Handleable interface */
    //@{
    /** Return a ref count pointer to a newly created object */
    virtual RCP<LinearSolverBase<double> > getRcp() 
    {return rcp(this);}
    //@}


  protected:

  private:
    
    void setupML(Epetra_RowMatrix* A) const ;

    /** Aztec options */
    mutable Array<int> options_;

    /** Aztec parameters */
    mutable Array<double> parameters_;

    /** Flag indicating whether we are using ML preconditioning */
    bool useML_;

    /** Flag indicating whether we are using Ifpack preconditioning */
    bool useIfpack_;

    /** Flag indicating whether we are using a user-defined preconditioner */
    bool useUserPrec_;

    /** Flag indicating whether we are doing a recursive solve
     *  with aztec (i.e., using recursiveIterate) */
    bool aztec_recursive_iterate_;

    /** Parameter list for preconditioner */
    mutable ParameterList precParams_;

    /** User-defined preconditioner object */
    mutable RCP<Epetra_Operator> userPrec_;

    /** Aztec status */
    mutable Array<double> aztec_status;

    /** Aztec proc_config */
    mutable Array<int> aztec_proc_config;


    /** Map from parameter name to AZTEC parameter identifier */
    static Teuchos::map<string,int>& paramMap() 
    {static Teuchos::map<string,int> rtn; return rtn;}
    
    /** Initialize the map from parameter names to AZTEC parameter ID codes */
    static void initParamMap();

    
  };
  
}

#endif
