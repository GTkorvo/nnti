
// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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

#ifndef ANASAZI_BLOCK_KRYLOV_SCHUR_SOLMGR_HPP
#define ANASAZI_BLOCK_KRYLOV_SCHUR_SOLMGR_HPP

/*! \file AnasaziBlockKrylovSchurSolMgr.hpp
 *  \brief The Anasazi::BlockKrylovSchurSolMgr provides a powerful solver manager for the BlockKrylovSchur eigensolver.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziModalSolverUtils.hpp"

#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestOrderedResNorm.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "Teuchos_BLAS.hpp"

/** \example BlockKrylovSchur/BlockKrylovSchurEpetraEx.cpp
    This is an example of how to use the Anasazi::BlockKrylovSchurSolMgr solver manager.
*/

/*! \class Anasazi::BlockKrylovSchurSolMgr
 *
 *  \brief The Anasazi::BlockKrylovSchurSolMgr provides a powerful and fully-featured solver manager over the BlockKrylovSchur eigensolver.
 *
 * This solver manager implements a hard-locking mechanism, whereby eigenpairs designated to be locked are moved from the eigensolver and placed in
 * auxiliary storage. The eigensolver is then restarted and continues to iterate, always orthogonal to the locked eigenvectors.

 \ingroup anasazi_solver_framework

 \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
 */

namespace Anasazi {

template<class ScalarType, class MV, class OP>
class BlockKrylovSchurSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MT;
    
  public:

  //! @name Constructors/Destructor
  //@{ 

  /*! \brief Basic constructor for BlockKrylovSchurSolMgr.
   *
   * This constructor accepts the Eigenproblem to be solved in addition
   * to a parameter list of options for the solver manager. These options include the following:
   *   - "Which" - a \c string specifying the desired eigenvalues: SM, LM, SR or LR. Default: "LM"
   *   - "Block Size" - a \c int specifying the block size to be used by the underlying block Krylov-Schur solver. Default: 1
   *   - "Num Blocks" - a \c int specifying the number of blocks allocated for the Krylov basis. Default: 3*nev
   *   - "Extra NEV Blocks" - a \c int specifying the number of extra blocks the solver should keep in addition to those
   *   -  required to compute the number of eigenvalues requested.  Default: 0
   *   - "Maximum Restarts" - a \c int specifying the maximum number of restarts the underlying solver is allowed to perform. Default: 20
   *   - "Orthogonalization" - a \c string specifying the desired orthogonalization:  DGKS and SVQB. Default: "SVQB"
   *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Anasazi::Errors
   *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: machine precision.
   *   - "Relative Convergence Tolerance" - a \c bool specifying whether residuals norms should be scaled by their eigenvalues for the purposing of deciding convergence. Default: true
   */
  BlockKrylovSchurSolMgr( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
                             Teuchos::ParameterList &pl );

  //! Destructor.
  virtual ~BlockKrylovSchurSolMgr() {};
  //@}
  
  //! @name Accessor methods
  //@{ 

  Eigenproblem<ScalarType,MV,OP>& getProblem() const {
    return *_problem;
  }

  /*! \brief Return the Ritz values from the most recent solve.
   */
  std::vector<Value<ScalarType> > getRitzValues() const {
    std::vector<Value<ScalarType> > ret( _ritzValues );
    return ret;
  }

  //@}

  //! @name Solver application methods
  //@{ 
    
  /*! \brief This method performs possibly repeated calls to the underlying eigensolver's iterate() routine
   * until the problem has been solved (as decided by the solver manager) or the solver manager decides to 
   * quit.
   *
   * This method calls BlockKrylovSchur::iterate(), which will return either because a specially constructed status test evaluates to ::Passed
   * or an exception is thrown.
   *
   * A return from BlockKrylovSchur::iterate() signifies one of the following scenarios:
   *    - the maximum number of restarts has been exceeded. In this scenario, the solver manager will place\n
   *      all converged eigenpairs into the eigenproblem and return ::Unconverged.
   *    - the locking conditions have been met. In this scenario, some of the current eigenpairs will be removed\n
   *      from the eigensolver and placed into auxiliary storage. The eigensolver will be restarted with the remaining part of the Krylov subspace\n
   *      and some random information to replace the removed subspace.
   *    - global convergence has been met. In this case, the most significant NEV eigenpairs in the solver and locked storage  \n
   *      have met the convergence criterion. (Here, NEV refers to the number of eigenpairs requested by the Eigenproblem.)    \n
   *      In this scenario, the solver manager will return ::Converged.
   *
   * \returns ::ReturnType specifying:
   * <ul>
   *     - ::Converged: the eigenproblem was solved to the specification required by the solver manager.
   *     - ::Unconverged: the eigenproblem was not solved to the specification desired by the solver manager.
   * </ul>
  */
  ReturnType solve();
  //@}

  private:
  Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > _problem;
  Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > _sort;

  string _whch, _ortho; 

  MagnitudeType _convtol;
  int _maxRestarts;
  bool _relconvtol,_conjSplit;
  int _blockSize, _numBlocks, _stepSize, _nevBlocks, _xtra_nevBlocks;
  int _verbosity;

  std::vector<Value<ScalarType> > _ritzValues;
};


// Constructor
template<class ScalarType, class MV, class OP>
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::BlockKrylovSchurSolMgr( 
        const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
        Teuchos::ParameterList &pl ) : 
  _problem(problem),
  _whch("LM"),
  _ortho("SVQB"),
  _convtol(0),
  _maxRestarts(20),
  _relconvtol(true),
  _conjSplit(false),
  _blockSize(0),
  _numBlocks(0),
  _stepSize(0),
  _nevBlocks(0),
  _xtra_nevBlocks(0),
  _verbosity(Anasazi::Errors)
{
  TEST_FOR_EXCEPTION(_problem == Teuchos::null,               std::invalid_argument, "Problem not given to solver manager.");
  TEST_FOR_EXCEPTION(!_problem->isProblemSet(),               std::invalid_argument, "Problem not set.");
  TEST_FOR_EXCEPTION(_problem->getInitVec() == Teuchos::null, std::invalid_argument, "Problem does not contain initial vectors to clone from.");

  const int nev = _problem->getNEV();

  // convergence tolerance
  _convtol = pl.get("Convergence Tolerance",MT::prec());
  _relconvtol = pl.get("Relative Convergence Tolerance",_relconvtol);
  
  // maximum number of restarts
  _maxRestarts = pl.get("Maximum Restarts",_maxRestarts);

  // block size: default is 1
  _blockSize = pl.get("Block Size",1);
  TEST_FOR_EXCEPTION(_blockSize <= 0, std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: \"Block Size\" must be strictly positive.");

  // set the number of blocks we need to save to compute the nev eigenvalues of interest.
  _xtra_nevBlocks = pl.get("Extra NEV Blocks",0);
  if (nev%_blockSize) {
    _nevBlocks = nev/_blockSize + _xtra_nevBlocks + 1;
  } else {
    _nevBlocks = nev/_blockSize + _xtra_nevBlocks;
  }

  _numBlocks = pl.get("Num Blocks",3*_nevBlocks);
  TEST_FOR_EXCEPTION(_numBlocks <= _nevBlocks, std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: \"Num Blocks\" must be strictly positive and large enough to compute the requested eigenvalues.");

  TEST_FOR_EXCEPTION(_numBlocks*_blockSize > MVT::GetVecLength(*_problem->getInitVec()),
                     std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: Potentially impossible orthogonality requests. Reduce basis size.");
  
  // step size: the default is _maxRestarts*_numBlocks, so that Ritz values are only computed every restart.
  if (_maxRestarts) {
    _stepSize = pl.get("Step Size", (_maxRestarts+1)*(_numBlocks+1));
  } else {
    _stepSize = pl.get("Step Size", _numBlocks+1);
  }
  TEST_FOR_EXCEPTION(_stepSize < 1, std::invalid_argument,
                     "Anasazi::BlockKrylovSchurSolMgr: \"Step Size\" must be strictly positive.");

  // get the sort manager
  if (pl.isParameter("Sort Manager")) {
    _sort = Teuchos::getParameter<Teuchos::RefCountPtr<Anasazi::SortManager<ScalarType,MV,OP> > >(pl,"Sort Manager");
  } else {
    // which values to solve for
    _whch = pl.get("Which",_whch);
    if (_whch != "SM" && _whch != "LM" && _whch != "SR" && _whch != "LR" && _whch != "SI" && _whch != "LI") {
      _whch = "LM";
    }
    _sort = Teuchos::rcp( new BasicSort<ScalarType,MV,OP>(_whch) );
  }

  // which orthogonalization to use
  _ortho = pl.get("Orthogonalization",_ortho);
  if (_ortho != "DGKS" && _ortho != "SVQB") {
    _ortho = "SVQB";
  }

  // verbosity level
  _verbosity = pl.get("Verbosity", _verbosity);

}


// solve()
template<class ScalarType, class MV, class OP>
ReturnType 
BlockKrylovSchurSolMgr<ScalarType,MV,OP>::solve() {

  const int nev = _problem->getNEV();
  ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
  ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();

  Teuchos::BLAS<int,ScalarType> blas;

  //////////////////////////////////////////////////////////////////////////////////////
  // Output manager
  Teuchos::RefCountPtr<BasicOutputManager<ScalarType> > printer = Teuchos::rcp( new BasicOutputManager<ScalarType>(_verbosity) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Status tests
  //
  // convergence
  Teuchos::RefCountPtr<StatusTestOrderedResNorm<ScalarType,MV,OP> > convtest 
    = Teuchos::rcp( new StatusTestOrderedResNorm<ScalarType,MV,OP>(_sort,_convtol,nev,StatusTestOrderedResNorm<ScalarType,MV,OP>::RITZRES_2NORM,_relconvtol) );

  // printing StatusTest
  Teuchos::RefCountPtr<StatusTestOutput<ScalarType,MV,OP> > outputtest
    = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer,convtest,1,Passed ) );

  //////////////////////////////////////////////////////////////////////////////////////
  // Orthomanager
  //
  Teuchos::RefCountPtr<OrthoManager<ScalarType,MV> > ortho; 
  if (_ortho=="SVQB") {
    ortho = Teuchos::rcp( new SVQBOrthoManager<ScalarType,MV,OP>(_problem->getM()) );
  } else if (_ortho=="DGKS") {
    ortho = Teuchos::rcp( new BasicOrthoManager<ScalarType,MV,OP>(_problem->getM()) );
  } else {
    TEST_FOR_EXCEPTION(_ortho!="SVQB"&&_ortho!="DGKS",std::logic_error,"Anasazi::BlockKrylovSchurSolMgr::solve(): Invalid orthogonalization type.");
  }
  
  // utils
  ModalSolverUtils<ScalarType,MV,OP> msutils(printer);

  //////////////////////////////////////////////////////////////////////////////////////
  // Parameter list
  Teuchos::ParameterList plist;
  plist.set("Block Size",_blockSize);
  plist.set("Num Blocks",_numBlocks);
  plist.set("Step Size",_stepSize);
  plist.set("Print Number of Ritz Values",_nevBlocks*_blockSize);

  //////////////////////////////////////////////////////////////////////////////////////
  // BlockKrylovSchur solver
  Teuchos::RefCountPtr<BlockKrylovSchur<ScalarType,MV,OP> > bks_solver 
    = Teuchos::rcp( new BlockKrylovSchur<ScalarType,MV,OP>(_problem,_sort,printer,outputtest,ortho,plist) );
  // set any auxiliary vectors defined in the problem
  Teuchos::RefCountPtr< const MV > probauxvecs = _problem->getAuxVecs();
  if (probauxvecs != Teuchos::null) {
    bks_solver->setAuxVecs( Teuchos::tuple< Teuchos::RefCountPtr<const MV> >(probauxvecs) );
  }

  // go ahead and initialize the solution to nothing in case we throw an exception
  Eigensolution<ScalarType,MV> sol;
  sol.numVecs = 0;
  _problem->setSolution(sol);

  int numRestarts = 0;
  int cur_nevBlocks = 0;

  // tell bks_solver to iterate
  while (1) {
    try {
      bks_solver->iterate();
  
      // check convergence first
      if (convtest->getStatus() == Passed ) {
        // we have convergence
        // convtest->whichVecs() tells us which vectors from solver state are the ones we want
        // convtest->howMany() will tell us how many
        break;
      }
      // check for restarting:  this is for the Hermitian case, or non-Hermitian conjugate split situation.
      // --> for the Hermitian case the current subspace dimension needs to match the maximum subspace dimension
      // --> for the non-Hermitican case:
      //     --> if a conjugate pair was detected in the previous restart then the current subspace dimension needs to match the
      //         maximum subspace dimension (the BKS solver keeps one extra vector if the problem is non-Hermitian).
      //     --> if a conjugate pair was not detected in the previous restart then the current subspace dimension will be one less
      //         than the maximum subspace dimension.
      else if ( (bks_solver->getCurSubspaceDim() == bks_solver->getMaxSubspaceDim()) ||
		(!_problem->isHermitian() && !_conjSplit && (bks_solver->getCurSubspaceDim()+1 == bks_solver->getMaxSubspaceDim())) ) {

        if ( numRestarts >= _maxRestarts ) {
          break; // break from while(1){bks_solver->iterate()}
        }
        numRestarts++;

        printer->stream(Debug) << " Performing restart number " << numRestarts << " of " << _maxRestarts << endl << endl;

	// Update the Schur form of the projected eigenproblem, then sort it.
	if (!bks_solver->isSchurCurrent())
	  bks_solver->computeSchurForm( true );
	
	// Get the most current Ritz values before we continue.
	_ritzValues = bks_solver->getRitzValues();

	// Get the state.
	BlockKrylovSchurState<ScalarType,MV> oldState = bks_solver->getState();
	
	// Get the current dimension of the factorization
	int curDim = oldState.curDim;

	// Determine if the storage for the nev eigenvalues of interest splits a complex conjugate pair.
	std::vector<int> ritzIndex = bks_solver->getRitzIndex();
	if (ritzIndex[_nevBlocks*_blockSize-1]==1) {
	  _conjSplit = true;
	  cur_nevBlocks = _nevBlocks*_blockSize+1;
	} else {
	  _conjSplit = false;
	  cur_nevBlocks = _nevBlocks*_blockSize;
	}
	
	// Update the Krylov-Schur decomposition

	// Get a view of the Schur vectors of interest.
	Teuchos::SerialDenseMatrix<int,ScalarType> Qnev(Teuchos::View, *(oldState.Q), curDim, cur_nevBlocks);

	// Get a view of the current Krylov basis.
	std::vector<int> curind( curDim );
	for (int i=0; i<curDim; i++) { curind[i] = i; }
	Teuchos::RefCountPtr<const MV> basistemp = MVT::CloneView( *(oldState.V), curind );

	// Create a vector for the new Krylov basis.
	// Need cur_nevBlocks for the updated factorization and another block for the current factorization residual block (F).
	Teuchos::RefCountPtr<MV> newV = MVT::Clone( *(oldState.V), cur_nevBlocks + _blockSize );
	curind.resize(cur_nevBlocks);
	Teuchos::RefCountPtr<MV> tmp_newV = MVT::CloneView( *newV, curind );
	
	// Compute the new Krylov basis.
	// --> Vnew = V*Qnev
	MVT::MvTimesMatAddMv( one, *basistemp, Qnev, zero, *tmp_newV );

	// Move the current factorization residual block (F) to the last block of newV.
	curind.resize(_blockSize);
	for (int i=0; i<_blockSize; i++) { curind[i] = curDim + i; }
        Teuchos::RefCountPtr<const MV> oldF = MVT::CloneView( *(oldState.V), curind );
	for (int i=0; i<_blockSize; i++) { curind[i] = cur_nevBlocks + i; }
	MVT::SetBlock( *oldF, curind, *newV );
	
	// Update the Krylov-Schur quasi-triangular matrix.
	
	// Create storage for the new Schur matrix of the Krylov-Schur factorization
	// Copy over the current quasi-triangular factorization of oldState.H which is stored in oldState.S.
	Teuchos::SerialDenseMatrix<int,ScalarType> oldS(Teuchos::View, *(oldState.S), cur_nevBlocks+_blockSize, cur_nevBlocks);
	Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > newH = 
	  Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( oldS ) );

	// Get a view of the B block of the current factorization
	Teuchos::SerialDenseMatrix<int,ScalarType> oldB(Teuchos::View, *(oldState.H), _blockSize, _blockSize, curDim, curDim-_blockSize);

	// Get a view of the a block row of the Schur vectors.
	Teuchos::SerialDenseMatrix<int,ScalarType> subQ(Teuchos::View, *(oldState.Q), _blockSize, cur_nevBlocks, curDim-_blockSize);
	
	// Get a view of the new B block of the updated Krylov-Schur factorization
	Teuchos::SerialDenseMatrix<int,ScalarType> newB(Teuchos::View, *newH,  _blockSize, cur_nevBlocks, cur_nevBlocks);

	// Compute the new B block.
	blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, _blockSize, cur_nevBlocks, _blockSize, one, 
		   oldB.values(), oldB.stride(), subQ.values(), subQ.stride(), zero, newB.values(), newB.stride() );
 
	// Set the new state and initialize the solver.
        BlockKrylovSchurState<ScalarType,MV> newstate;
        newstate.V = newV;
	newstate.H = newH;
	newstate.curDim = cur_nevBlocks;
        bks_solver->initialize(newstate);

      }
      else {
        TEST_FOR_EXCEPTION(true,std::logic_error,"Anasazi::BlockKrylovSchurSolMgr::solve(): Invalid return from bks_solver::iterate().");
      }
    }
    catch (std::exception e) {
      printer->stream(Errors) << "Error! Caught exception in BlockKrylovSchur::iterate() at iteration " << bks_solver->getNumIters() << endl 
                              << e.what() << endl;
      throw;
    }
  }

  // Get the most current Ritz values before we return
  _ritzValues = bks_solver->getRitzValues();
  
  sol.numVecs = convtest->howMany();
  if (sol.numVecs > 0) {
    sol.index = bks_solver->getRitzIndex();
    sol.Evals = bks_solver->getRitzValues();
    // Check to see if conjugate pair is on the boundary.
    if (sol.index[sol.numVecs-1]==1) {
      sol.numVecs++;
      sol.Evals.resize(sol.numVecs);
      sol.index.resize(sol.numVecs);
      bks_solver->setNumRitzVectors(sol.numVecs);
    } else {
      sol.Evals.resize(sol.numVecs);
      sol.index.resize(sol.numVecs);
      bks_solver->setNumRitzVectors(sol.numVecs);
    }
    bks_solver->computeRitzVectors();
    sol.Evecs = MVT::CloneCopy( *(bks_solver->getRitzVectors()) );
    sol.Espace = sol.Evecs;
  } 

  // print final summary
  bks_solver->currentStatus(printer->stream(FinalSummary));

  // print timing information
  Teuchos::TimeMonitor::summarize(printer->stream(TimingDetails));

  _problem->setSolution(sol);
  printer->stream(Debug) << "Returning " << sol.numVecs << " eigenpairs to eigenproblem." << endl;

  if (sol.numVecs < nev) {
    return Unconverged; // return from BlockKrylovSchurSolMgr::solve() 
  }
  return Converged; // return from BlockKrylovSchurSolMgr::solve() 
}


} // end Anasazi namespace

#endif /* ANASAZI_BLOCK_KRYLOV_SCHUR_SOLMGR_HPP */
