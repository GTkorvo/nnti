// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


/** 
    \class BelosKrylov
    \brief Provides interface for using ROL::Vector with Belos solvers
                
    \author Created by Greg von Winckel
*/


#ifndef ROL_BELOS_KRYLOV_HPP
#define ROL_BELOS_KRYLOV_HPP

#include "BelosSolverFactory.hpp"   
#include "BelosSolverManager.hpp"   

#include "ROL_Krylov.hpp"
#include "ROL_BelosVector.hpp"
#include "ROL_BelosOperator.hpp"

namespace ROL {

    template<class Real>
    class BelosKrylov : public Krylov<Real> {
 
        typedef Real               ST;
        typedef LinearOperator<ST> OP;
        typedef Vector<ST>         MV;

        // For testing
	typedef Belos::MultiVecTraits<ST,MV>    MVT;
        typedef Belos::OperatorTraits<ST,MV,OP> OPT;

        private:

            Belos::SolverFactory<ST,MV,OP> factory_;
            Teuchos::RCP<Belos::SolverManager<ST,MV,OP>> solver_;
            Teuchos::RCP<Belos::LinearProblem<ST,MV,OP>> problem_;  

        public:
           
            /// \brief Create a Belos solver 
            BelosKrylov(Teuchos::ParameterList &parlist) : 
                problem_(Teuchos::rcp(new Belos::LinearProblem<ST,MV,OP>)) {

                auto solverParams = Teuchos::rcp(new Teuchos::ParameterList());

                // Options likely to be of interest include CG, MINRES, GMRES, GCRODR, and RCG
                auto blockSize          = 1; // Only support single solution & single RHS for now 
                auto solverName         = parlist.get("Belos Krylov Method","MINRES");  
                auto maxit              = parlist.get("Maximum Number of Krylov Iterations",50);
                auto abstol             = parlist.get("Absolute Krylov Tolerance",1.e-4);
                auto numVectors         = parlist.get("Number of Stored Vectors",3);
                
                solverParams->setName("Belos input parameters"); 
                solverParams->set("Block Size",blockSize);
                solverParams->set("Maximum Iterations",maxit);
                solverParams->set("Convergence Tolerance",abstol);  
                solverParams->set("Num Blocks",numVectors);

                solver_ = factory_.create(solverName,solverParams);                 
            }


            /// \brief Compute solution vector
            void run( MV &x, OP& A, const MV &b, OP &M, int &iter, int &flag )  {

                // Need to get RCPs for x,A,b, and M
                Teuchos::RCP<MV>       xp = Teuchos::rcpFromRef(x);
                Teuchos::RCP<OP>       Ap = Teuchos::rcpFromRef(A);
                Teuchos::RCP<const MV> bp = Teuchos::rcpFromRef(b);
                Teuchos::RCP<OP>       Mp = Teuchos::rcpFromRef(M);

                problem_->setOperator(Ap);
                problem_->setLeftPrec(Mp);
                problem_->setProblem(xp,bp);

                solver_->setProblem(problem_);

                flag = static_cast<int>(solver_->solve());

                iter = solver_->getNumIters();
            }
    };
}

#endif
