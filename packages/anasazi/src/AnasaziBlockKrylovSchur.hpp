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

/*! \file AnasaziBlockKrylovSchur.hpp
  \brief Implementation of a block Krylov-Schur eigensolver.
*/

// TODO: the documentation here needs to be made rigorous
// in particular, getState() and initialize() need to exactly describe their 
// input and output

#ifndef ANASAZI_BLOCK_KRYLOV_SCHUR_HPP
#define ANASAZI_BLOCK_KRYLOV_SCHUR_HPP

#include "AnasaziTypes.hpp"

#include "AnasaziEigensolver.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "AnasaziOrthoManager.hpp"
#include "AnasaziModalSolverUtils.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

#if defined(HAVE_COMPLEX)
#define ANSZI_CPLX_CLASS std::complex
#elif  defined(HAVE_COMPLEX_H)
#define ANSZI_CPLX_CLASS ::complex
#endif

/*!     \class Anasazi::BlockKrylovSchur
  
  \brief This class implements the block Krylov-Schur iteration,
  for solving eigenvalue problems.
  
  This method is a block version of the iteration presented by G.W. Stewart 
  in "A Krylov-Schur Algorithm for Large Eigenproblems", 
  SIAM J. Matrix Anal. Appl., Vol 23(2001), No. 3, pp. 601-614.
  
  \ingroup anasazi_solvers

  \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {

  //! @name BlockKrylovSchur Structures 
  //@{ 

  /** \brief Structure to contain pointers to BlockKrylovSchur state variables.
   *
   * This struct is utilized by BlockKrylovSchur::initialize() and BlockKrylovSchur::getState().
   */
  template <class ScalarType, class MulVec>
  struct BlockKrylovSchurState {
    int curDim;
    Teuchos::RefCountPtr<const MulVec> V;
    Teuchos::RefCountPtr<const MulVec> X;
    Teuchos::RefCountPtr<const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > T;
    Teuchos::RefCountPtr<const Teuchos::SerialDenseMatrix<int,ScalarType> > H;
    BlockKrylovSchurState() : curDim(0), V(Teuchos::null),
                              X(Teuchos::null), T(Teuchos::null), 
			      H(Teuchos::null) {}
  };

  //@}

  //! @name BlockKrylovSchur Exceptions
  //@{ 

  /** \brief BlockKrylovSchurInitFailure is thrown when the BlockKrylovSchur solver is unable to
   * generate an initial iterate in the BlockKrylovSchur::initialize() routine. 
   *
   * This exception is thrown from the BlockKrylovSchur::initialize() method, which is
   * called by the user or from the BlockKrylovSchur::iterate() method if isInitialized()
   * == \c false.
   *
   * In the case that this exception is thrown, 
   * BlockKrylovSchur::isInitialized() will be \c false and the user will need to provide
   * a new initial iterate to the solver.
   *
   */
  class BlockKrylovSchurInitFailure : public AnasaziError {public:
    BlockKrylovSchurInitFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};

  /** \brief BlockKrylovSchurOrthoFailure is thrown when the orthogonalization manager is
   * unable to orthogonalize the preconditioned residual against (a.k.a. \c H)
   * the current basis (a.k.a. \c V).
   *
   * This exception is thrown from the BlockKrylovSchur::iterate() method.
   *
   */
  class BlockKrylovSchurOrthoFailure : public AnasaziError {public:
    BlockKrylovSchurOrthoFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};
  
  //@}


  template <class ScalarType, class MV, class OP>
  class BlockKrylovSchur : public Eigensolver<ScalarType,MV,OP> { 
  public:
    //! @name Constructor/Destructor
    //@{ 
    
    //! %Anasazi::BlockKrylovSchur constructor.
    BlockKrylovSchur( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
                   const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sorter,
                   const Teuchos::RefCountPtr<OutputManager<ScalarType> > &printer,
                   const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > &tester,
                   const Teuchos::RefCountPtr<OrthoManager<ScalarType,MV> > &ortho,
                   Teuchos::ParameterList &params 
                 );
    
    //! %Anasazi::BlockKrylovSchur destructor.
    virtual ~BlockKrylovSchur() {};
    //@}


    //! @name Solver methods
    //@{ 
    
    /*! \brief This method performs %BlockKrylovSchur iterations until the status
     * test indicates the need to stop or an error occurs (in which case, an
     * exception is thrown).
     *
     * iterate() will first determine whether the solver is inintialized; if
     * not, it will call initialize() using default arguments. After
     * initialization, the solver performs %BlockKrylovSchur iterations until the
     * status test evaluates as Passed, at which point the method returns to
     * the caller. 
     *
     */
    void iterate();

    /*! \brief Initialize the solver to an iterate, optionally providing the
     * Ritz values, residual, and search direction.
     *
     * The %BlockKrylovSchur eigensolver contains a certain amount of state
     * relating to the current eigenvectors, including the current residual,
     * the current Krylov basis, and the images of these under the operators.
     *
     * initialize() gives the user the opportunity to manually set these,
     * although this must be done with caution, abiding by the rules given
     * below. All notions of orthogonality and orthonormality are derived from
     * the inner product specified by the orthogonalization manager.
     *
     * \post 
     * <li>isInitialized() == true (see post-conditions of isInitialize())
     *
     * The user has the option of specifying any component of the state using
     * initialize(). However, these arguments are assumed to match the
     * post-conditions specified under isInitialized(). Any component of the
     * state (i.e., KX) not given to initialize() will be generated.
     *
     */
    void initialize(BlockKrylovSchurState<ScalarType,MV> state);
    void initialize();

    /*! \brief Indicates whether the solver has been initialized or not.
     *
     * \return bool indicating the state of the solver.
     * \post
     * <ul>
     * <li> finish
     * </ul>
     */
    bool isInitialized() { return _initialized; }

    /*! \brief Get the current state of the eigensolver.
     * 
     * The data is only valid if isInitialized() == \c true. 
     *
     * The data for the preconditioned residual is only meaningful in the
     * scenario that the solver throws an ::BlockKrylovSchurRitzFailure exception
     * during iterate().
     *
     * \returns A BlockKrylovSchurState object containing const pointers to the current
     * solver state.
     */
    BlockKrylovSchurState<ScalarType,MV> getState() const {
      BlockKrylovSchurState<ScalarType,MV> state;
      state.curDim = _curDim;
      state.V = _V;
      state.X = _X;
      state.H = _H;
      state.T = Teuchos::rcp(new std::vector<MagnitudeType>(_ritzvalues));
      return state;
    }
    
    //@}


    //! @name Status methods
    //@{ 

    //! \brief Get the current iteration count.
    int getNumIters() const { return(_iter); };

    //! \brief Reset the iteration count.
    void resetNumIters() { _iter=0; };

    //! \brief Get the current approximate eigenvectors.
    Teuchos::RefCountPtr<const MV> getEvecs() {return _X;}

    //! \brief Get the residual vectors.
    Teuchos::RefCountPtr<const MV> getResidualVecs() {return Teuchos::null;}

    /*! \brief Get the current eigenvalue estimates.
     *
     *  \return A vector of length getBlockSize() containing the eigenvalue
     *  estimates associated with the current iterate.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getEigenvalues() { 
      std::vector<MagnitudeType> ret = _ritzvalues;
      ret.resize(_blockSize);
      return ret;
    }

    /*! \brief Get the Ritz values for the previous iteration.
     *
     *  \return A vector of length not exceeding the maximum dimension of the subspace 
     *  containing the Ritz values from the previous projected eigensolve.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRitzValues() { 
      std::vector<MagnitudeType> ret = _ritzvalues;
      ret.resize(_curDim);
      return ret;
    }

    /*! \brief Get the current residual norms
     *
     *  \return A vector of length blockSize containing the norms of the
     *  residuals, with respect to the orthogonalization manager norm() method.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getResNorms()    {return _Rnorms;}


    /*! \brief Get the current residual 2-norms
     *
     *  \return A vector of length blockSize containing the 2-norms of the
     *  residuals. 
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRes2Norms()   {return _R2norms;}

    /*! \brief Get the current ritz residual 2-norms
     *
     *  \return A vector of length blockSize containing the 2-norms of the
     *  ritz residuals.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRitzRes2Norms() {return _ritzresiduals;}

    //! @name Accessor routines
    //@{ 


    //! Get a constant reference to the eigenvalue problem.
    const Eigenproblem<ScalarType,MV,OP>& getProblem() const { return(*_problem); };


    /*! \brief Set the blocksize and number of blocks to be used by the
     * iterative solver in solving this eigenproblem.
     *  
     *  Changing either the block size or the number of blocks will reset the
     *  solver to an uninitialized state.
     */
    void setSize(int blockSize, int numBlocks);

    //! \brief Set the blocksize. 
    void setBlockSize(int blockSize);

    //! \brief Set the step size. 
    void setStepSize(int stepSize);

    //! \brief Get the step size. 
    int getStepSize() {return _stepSize;}

    //! Get the blocksize to be used by the iterative solver in solving this eigenproblem.
    int getBlockSize() const { return(_blockSize); }

    /*! \brief Get the dimension of the search subspace used to generate the current eigenvectors and eigenvalues.
     *
     *  \return An integer specifying the rank of the Krylov subspace currently in use by the eigensolver. If isInitialized() == \c false, 
     *  the return is 0. Otherwise, it will be some strictly positive multiple of getBlockSize().
     */
    int getCurSubspaceDim() {
      if (!_initialized) return 0;
      return _curDim;
    }

    //! Get the maximum dimension allocated for the search subspace. For %BlockKrylovSchur, this always returns numBlocks*blockSize.
    int getMaxSubspaceDim() {return _blockSize*_numBlocks;}


    /*! \brief Set the auxiliary vectors for the solver.
     *
     *  Because the current iterate X and search direction P cannot be assumed
     *  orthogonal to the new auxiliary vectors, a call to setAuxVecs() may
     *  reset the solver to the uninitialized state. This happens only in the
     *  case where the user requests full orthogonalization when the solver is
     *  in an initialized state with full orthogonalization disabled.
     *
     *  In order to preserve the current iterate, the user will need to extract
     *  it from the solver using getEvecs(), orthogonalize it against the new
     *  auxiliary vectors, and manually reinitialize the solver using
     *  initialize().
     */
    void setAuxVecs(const Teuchos::Array<Teuchos::RefCountPtr<const MV> > &auxvecs);

    //! Get the auxiliary vectors for the solver.
    Teuchos::Array<Teuchos::RefCountPtr<const MV> > getAuxVecs() const {return _auxVecs;}

    //@}

    //! @name Output methods
    //@{ 
    
    //! This method requests that the solver print out its current status to screen.
    void currentStatus(ostream &os);

    //@}

  private:
    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    typedef typename std::vector<ScalarType>::iterator STiter;
    typedef typename std::vector<MagnitudeType>::iterator MTiter;
    const MagnitudeType ONE;  
    const MagnitudeType ZERO; 
    const MagnitudeType NANVAL;
    //
    // Internal structs
    //
    struct CheckList {
      bool checkV;
      bool checkX;
      bool checkAux;
      CheckList() : checkV(false), checkX(false), checkAux(false) {};
    };
    //
    // Internal methods
    //
    string accuracyCheck(const CheckList &chk, const string &where) const;
    //
    // Classes inputed through constructor that define the eigenproblem to be solved.
    //
    const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> >     _problem;
    const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> >      _sm;
    const Teuchos::RefCountPtr<OutputManager<ScalarType> >          _om;
    const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> >       _tester;
    const Teuchos::RefCountPtr<OrthoManager<ScalarType,MV> >        _orthman;
    //
    // Information obtained from the eigenproblem
    //
    Teuchos::RefCountPtr<OP> _Op;
    //
    // Internal utilities class required by eigensolver.
    //
    ModalSolverUtils<ScalarType,MV,OP> _MSUtils;
    //
    // Internal timers
    //
    Teuchos::RefCountPtr<Teuchos::Time> _timerOp, _timerSortEval,
                                        _timerCompSF, _timerSortSF,
                                        _timerCompEvec, _timerQRFact,
                                        _timerOrtho, _timerTotal;
    //
    // Counters
    //
    int _count_ApplyOp;

    //
    // Algorithmic parameters.
    //
    // _blockSize is the solver block size; it controls the number of eigenvectors that 
    // we compute, the number of residual vectors that we compute, and therefore the number
    // of vectors added to the basis on each iteration.
    int _blockSize;
    // _numBlocks is the size of the allocated space for the Krylov basis, in blocks.
    int _numBlocks; 
    // _stepSize dictates how many iterations are performed before eigenvectors and eigenvalues
    // are computed again
    int _stepSize;
    
    // 
    // Current solver state
    //
    // _initialized specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of _initialized, please see documentation for initialize()
    bool _initialized;
    //
    // _isVecsCurrent describes whether the current eigenvectors correspond to the current basis
    bool _isVecsCurrent;
    // _curDim reflects how much of the current basis is valid 
    // NOTE: 0 <= _curDim <= _blockSize*_numBlocks
    // this also tells us how many of the values in _theta are valid Ritz values
    int _curDim;
    //
    // State Multivecs
    Teuchos::RefCountPtr<MV> _X, _V;
    //
    // Projected matrices
    //
    Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > _H;
    // 
    // Auxiliary vectors
    Teuchos::Array<Teuchos::RefCountPtr<const MV> > _auxVecs;
    int _numAuxVecs;
    //
    // Number of iterations that have been performed.
    int _iter;
    // 
    // Current eigenvalues, residual norms
    std::vector<MagnitudeType> _ritzvalues, _Rnorms, _R2norms, _ritzresiduals;
    //
    // Current Schur Error || Op*V - VR ||.
    MagnitudeType _schurError;

  };


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor
  template <class ScalarType, class MV, class OP>
  BlockKrylovSchur<ScalarType,MV,OP>::BlockKrylovSchur(
        const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
        const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sorter,
        const Teuchos::RefCountPtr<OutputManager<ScalarType> > &printer,
        const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > &tester,
        const Teuchos::RefCountPtr<OrthoManager<ScalarType,MV> > &ortho,
        Teuchos::ParameterList &params
        ) :
    ONE(Teuchos::ScalarTraits<MagnitudeType>::one()),
    ZERO(Teuchos::ScalarTraits<MagnitudeType>::zero()),
    NANVAL(Teuchos::ScalarTraits<MagnitudeType>::nan()),
    // problem, tools
    _problem(problem), 
    _sm(sorter),
    _om(printer),
    _tester(tester),
    _orthman(ortho),
    _Op(_problem->getOperator()),
    _MSUtils(_om),
    // timers, counters
    _timerOp(Teuchos::TimeMonitor::getNewTimer("Operation Op*x")),
    _timerSortEval(Teuchos::TimeMonitor::getNewTimer("Sorting eigenvalues")),
    _timerCompSF(Teuchos::TimeMonitor::getNewTimer("Computing Schur form")),
    _timerSortSF(Teuchos::TimeMonitor::getNewTimer("Sorting Schur form")),
    _timerCompEvec(Teuchos::TimeMonitor::getNewTimer("Computing eigenvectors")),
    _timerQRFact(Teuchos::TimeMonitor::getNewTimer("QR factorization")),
    _timerOrtho(Teuchos::TimeMonitor::getNewTimer("Orthogonalization")),
    _timerTotal(Teuchos::TimeMonitor::getNewTimer("Total time")),
    _count_ApplyOp(0),
    // internal data
    _blockSize(0),
    _numBlocks(0),
    _stepSize(0),
    _initialized(false),
    _curDim(0),
    _auxVecs( Teuchos::Array<Teuchos::RefCountPtr<const MV> >(0) ), 
    _numAuxVecs(0),
    _iter(0),
    _schurError(ONE)
  {     
    TEST_FOR_EXCEPTION(_problem == Teuchos::null,std::logic_error,
                       "Anasazi::BlockKrylovSchur::constructor: user specified null problem pointer.");
    TEST_FOR_EXCEPTION(_problem->isProblemSet() == false, std::logic_error,
                       "Anasazi::BlockKrylovSchur::constructor: user specified problem is not set.");

    // set the block size and allocate data
    int bs = params.get("Block Size", _problem->getNEV());
    int nb = params.get("Num Blocks", 1);
    setSize(bs,nb);
    int ss = params.get("Step Size",_numBlocks);
    setStepSize(ss);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size
  // This simply calls setSize(), modifying the block size while retaining the number of blocks.
  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::setBlockSize (int blockSize) 
  {
    setSize(blockSize,_numBlocks);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the step size.
  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::setStepSize (int stepSize)
  {
    TEST_FOR_EXCEPTION(stepSize <= 0, std::logic_error, "Anasazi::BlockKrylovSchur::setStepSize(): new step size must be positive and non-zero.");
    _stepSize = stepSize;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::setSize (int blockSize, int numBlocks) 
  {
    // This routine only allocates space; it doesn't not perform any computation
    // any change in size will invalidate the state of the solver.

    TEST_FOR_EXCEPTION(numBlocks <= 0 || blockSize <= 0, std::logic_error, "Anasazi::BlockKrylovSchur::setSize was passed a non-positive argument.");
    if (blockSize == _blockSize && numBlocks == _numBlocks) {
      // do nothing
      return;
    }

    _blockSize = blockSize;
    _numBlocks = numBlocks;

    Teuchos::RefCountPtr<const MV> tmp;
    // grab some Multivector to Clone
    // in practice, getInitVec() should always provide this, but it is possible to use a 
    // Eigenproblem with nothing in getInitVec() by manually initializing with initialize(); 
    // in case of that strange scenario, we will try to Clone from _V first, then resort to getInitVec()
    if (_V != Teuchos::null) { // this is equivalent to _blockSize > 0
      tmp = _V;
    }
    else {
      tmp = _problem->getInitVec();
      TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::logic_error,
                         "Anasazi::BlockKrylovSchur::setSize(): Eigenproblem did not specify initial vectors to clone from");
    }

    //////////////////////////////////
    // blockSize dependent
    //
    // grow/allocate vectors
    _Rnorms.resize(_blockSize,NANVAL);
    _R2norms.resize(_blockSize,NANVAL);
    //
    // clone multivectors off of tmp
    _X = MVT::Clone(*tmp,_blockSize);

    //////////////////////////////////
    // blockSize*numBlocks dependent
    //
    int newsd = _blockSize*_numBlocks;
    _ritzvalues.resize(newsd,NANVAL);
    _ritzresiduals.resize(newsd,NANVAL);
    _V = MVT::Clone(*tmp,newsd);
    _H = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newsd+_blockSize,newsd) );

    _initialized = false;
    _curDim = 0;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the auxiliary vectors
  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::setAuxVecs(const Teuchos::Array<Teuchos::RefCountPtr<const MV> > &auxvecs) {
    typedef typename Teuchos::Array<Teuchos::RefCountPtr<const MV> >::iterator tarcpmv;

    // set new auxiliary vectors
    _auxVecs = auxvecs;
    
    if (_om->isVerbosity( Debug ) ) {
      // Check almost everything here
      CheckList chk;
      chk.checkAux = true;
      _om->print( Debug, accuracyCheck(chk, ": in setAuxVecs()") );
    }

    _numAuxVecs = 0;
    for (tarcpmv i=_auxVecs.begin(); i != _auxVecs.end(); i++) {
      _numAuxVecs += MVT::GetNumberVecs(**i);
    }
    
    // If the solver has been initialized, X and P are not necessarily orthogonal to new auxiliary vectors
    if (_numAuxVecs > 0 && _initialized) {
      _initialized = false;
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  /* Initialize the state of the solver
   * 
   * POST-CONDITIONS:
   *
   * _V is orthonormal, orthogonal to _auxVecs, for first _curDim vectors
   * _ritzvalues contains Ritz w.r.t. V
   * X is Ritz vectors w.r.t. V
   */
  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::initialize(BlockKrylovSchurState<ScalarType,MV> state)
  {
    // NOTE: memory has been allocated by setBlockSize(). Use SetBlock below; do not Clone

    std::vector<int> bsind(_blockSize);
    for (int i=0; i<_blockSize; i++) bsind[i] = i;

    Teuchos::BLAS<int,ScalarType> blas;

    // in BlockKrylovSchur, V is primary
    // the order of dependence follows like so.
    // --init->               V,H
    //
    // inconsitent multivectors widths and lengths will not be tolerated, and
    // will be treated with exceptions.
    //
    std::string errstr("Anasazi::BlockKrylovSchur::initialize(): multivectors must have a consistent length and width.");

    // set up V,H: if the user doesn't specify these, ignore the rest
    if (state.V != Teuchos::null && state.H != Teuchos::null) {
      TEST_FOR_EXCEPTION( MVT::GetVecLength(*state.V) != MVT::GetVecLength(*_V),
                          std::logic_error, errstr );
      TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*state.V) < _blockSize,
                          std::logic_error, errstr );
      TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*state.V) > _blockSize*_numBlocks,
                          std::logic_error, errstr );

      _curDim = MVT::GetNumberVecs(*state.V);
      // pick an integral amount
      //_curDim = (int)(_curDim / _blockSize)*_blockSize;

      // this test is equivalent to _curDim == 0, but makes for more helpful output
      TEST_FOR_EXCEPTION( _curDim < _blockSize, BlockKrylovSchurInitFailure, "Anasazi::BlockKrylovSchur::initialize(): user-specified V must be >= blockSize.");
      // check size of H
      TEST_FOR_EXCEPTION(state.H->numRows() < _curDim || state.H->numCols() < _curDim, std::logic_error, errstr);

      // create index vector to copy over current basis vectors
      std::vector<int> nevind(_curDim);
      for (int i=0; i<_curDim; i++) nevind[i] = i;

      // copy basis vectors from the state into V
      MVT::SetBlock(*state.V,nevind,*_V);

      // get local view of V: view of first _curDim vectors
      Teuchos::RefCountPtr<MV> lclV;

      // generate lclV in case we need it below
      lclV = MVT::CloneView(*_V,nevind);

      // put data into _H
      Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > lclH;
      lclH = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*_H,_curDim,_curDim) );
      lclH->assign(*state.H);

      // X,ritzvalues require Ritz analisys; if we have to generate one of these, we might as well generate both
      if (state.X != Teuchos::null && state.T != Teuchos::null) {
        TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.X) != _blockSize,
                            std::logic_error, errstr );
        TEST_FOR_EXCEPTION(MVT::GetVecLength(*state.X) != MVT::GetVecLength(*_X),
                            std::logic_error, errstr );
        TEST_FOR_EXCEPTION((signed int)(state.T->size()) != _curDim,
                            std::logic_error, errstr );
        MVT::SetBlock(*state.X,bsind,*_X);
        std::copy(state.T->begin(),state.T->end(),_ritzvalues.begin());
      }
      else {
        // compute ritz vecs/vals
        Teuchos::SerialDenseMatrix<int,ScalarType> S(_curDim,_curDim);
        {
          Teuchos::TimeMonitor lcltimer( *_timerCompSF );
          int rank = _curDim;
          _MSUtils.directSolver(_curDim, *lclH, 0, &S, &_ritzvalues, &rank, 10);
          // we want all ritz values back
          TEST_FOR_EXCEPTION(rank != _curDim,BlockKrylovSchurInitFailure,
                             "Anasazi::BlockKrylovSchur::initialize(): Not enough Ritz vectors to initialize algorithm.");
        }
        // sort ritz pairs
        {
          Teuchos::TimeMonitor lcltimer( *_timerSortEval );

          std::vector<int> _order(_curDim);
          // make a ScalarType copy of ritzvalues
          std::vector<ScalarType> _ritzvalues_st(_curDim);
          std::copy(_ritzvalues.begin(),_ritzvalues.begin()+_curDim,_ritzvalues_st.begin());
          // sort it
          _sm->sort( NULL, _curDim, &(_ritzvalues_st[0]), &_order );   // don't catch exception
          // Put the sorted ritz values back into _ritzvalues
          for (int i=0; i<_curDim; i++) {
            _ritzvalues[i] = SCT::real(_ritzvalues_st[i]);
          }
          // Sort the primitive ritz vectors
          Teuchos::SerialDenseMatrix<int,ScalarType> copyS( S );
          for (int i=0; i<_curDim; i++) {
            blas.COPY(_curDim, copyS[_order[i]], 1, S[i], 1);
          }
        }
        Teuchos::SerialDenseMatrix<int,ScalarType> S1(Teuchos::View,S,_curDim,_blockSize);
        {
          //Teuchos::TimeMonitor lcltimer( *_timerLocal );

          // X <- lclV*S
          MVT::MvTimesMatAddMv( ONE, *lclV, S1, ZERO, *_X );
        }
      }

      // done with local pointers
      lclV = Teuchos::null;
      lclH = Teuchos::null;

      // Update the residual norms
      //      _orthman->norm(*_R,&_Rnorms);
      // Update the residual 2-norms 
      //      MVT::MvNorm(*_R,&_R2norms);

      _initialized = true;

      if (_om->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkV = true;
        chk.checkX = true;
        chk.checkAux = true;
        _om->print( Debug, accuracyCheck(chk, ": after initialize()") );
      }

      // Print information on current status
      if (_om->isVerbosity(Debug)) {
        currentStatus( _om->stream(Debug) );
      }
      else if (_om->isVerbosity(IterationDetails)) {
        currentStatus( _om->stream(IterationDetails) );
      }
    }
    else {
      // user did not specify a basis V
      // get vectors from problem or generate something, projectAndNormalize, call initialize() recursively
      Teuchos::RefCountPtr<const MV> ivec = _problem->getInitVec();
      TEST_FOR_EXCEPTION(ivec == Teuchos::null,std::logic_error,
                         "Anasazi::BlockDavdison::initialize(): Eigenproblem did not specify initial vectors to clone from");

      int lclDim = MVT::GetNumberVecs(*ivec);
      // pick the largest multiple of _blockSize
      lclDim = (int)(lclDim / _blockSize)*_blockSize;
      bool userand = false;
      if (lclDim < _blockSize) {
        // we need at least _blockSize vectors
        // use a random multivec
        userand = true;
        lclDim = _blockSize;
      }

      // make an index
      std::vector<int> dimind(lclDim);
      for (int i=0; i<lclDim; i++) dimind[i] = i;

      // alloc newV, newKV, newMV
      Teuchos::RefCountPtr<MV> newMV, 
                               newKV = MVT::Clone(*ivec,lclDim),
                               newV  = MVT::Clone(*ivec,lclDim);
      if (userand) {
        MVT::MvRandom(*newV);
      }
      else {
        // assign ivec to first part of newV
        MVT::SetBlock(*ivec,dimind,*newV);
      }

      // remove auxVecs from newV and normalize newV
      if (_auxVecs.size() > 0) {
        Teuchos::TimeMonitor lcltimer( *_timerOrtho );

        Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummy;
        int rank = _orthman->projectAndNormalize(*newV,dummy,Teuchos::null,_auxVecs);
        TEST_FOR_EXCEPTION(rank != lclDim,BlockKrylovSchurInitFailure,
                           "Anasazi::BlockKrylovSchur::initialize(): Couldn't generate initial basis of full rank.");
      }
      else {
        Teuchos::TimeMonitor lcltimer( *_timerOrtho );

        int rank = _orthman->normalize(*newV,Teuchos::null);
        TEST_FOR_EXCEPTION(rank != lclDim,BlockKrylovSchurInitFailure,
                           "Anasazi::BlockKrylovSchur::initialize(): Couldn't generate initial basis of full rank.");
      }

      // compute newKV
      {
        Teuchos::TimeMonitor lcltimer( *_timerOp );
        OPT::Apply(*_Op,*newV,*newKV);
        _count_ApplyOp += lclDim;
      }

      // generate KK
      Teuchos::RefCountPtr< Teuchos::SerialDenseMatrix<int,ScalarType> > H;
      H = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(lclDim,lclDim) );

      MVT::MvTransMv(ONE,*newV,*newKV,*H);

      // clear newHV,newMV
      newKV = newMV = Teuchos::null;

      // call myself recursively
      BlockKrylovSchurState<ScalarType,MV> newstate;
      newstate.V = newV;
      newstate.H = H;
      initialize(newstate);
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // initialize the solver with default state
  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::initialize()
  {
    BlockKrylovSchurState<ScalarType,MV> empty;
    initialize(empty);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Perform BlockKrylovSchur iterations until the StatusTest tells us to stop.
  template <class ScalarType, class MV, class OP>
  void BlockKrylovSchur<ScalarType,MV,OP>::iterate ()
  {
    //
    // Allocate/initialize data structures
    //
    if (_initialized == false) {
      initialize();
    }

    // as a data member, this would be redundant and require synchronization with
    // _blockSize and _numBlocks; we'll use a constant here.
    const int searchDim = _blockSize*_numBlocks;
    // we use this often enough...
    std::vector<int> bsind(_blockSize);
    for (int i=0; i<_blockSize; i++) { bsind[i] = i; }

    Teuchos::BLAS<int,ScalarType> blas;

    //
    // The projected matrices are part of the state, but the eigenvectors can be defined
    // locally.
    //    S = Local eigenvectors         (size: searchDim * searchDim
    Teuchos::SerialDenseMatrix<int,ScalarType> S( searchDim, searchDim );


    ////////////////////////////////////////////////////////////////
    // iterate until the status test tells us to stop.
    // also break if our basis is full
    while (_tester->checkStatus(this) != Passed && _curDim < searchDim) {

      _iter++;

      // Get the current part of the basis.
      std::vector<int> curind(_blockSize);
      for (int i=0; i<_blockSize; i++) curind[i] = _curDim + i;
      Teuchos::RefCountPtr<MV> Vnext = MVT::CloneView(*_V,curind);

      // Get a view of the previous vectors
      // this is used for orthogonalization and for computing V^H K H
      std::vector<int> prevind(_curDim);
      for (int i=0; i<_curDim; i++) prevind[i] = i;
      Teuchos::RefCountPtr<MV> Vprev = MVT::CloneView(*_V,prevind);

      // V has been extended, and H has been extended. Update basis dim and release all pointers.
      _curDim += _blockSize;
      Vnext = Teuchos::null;
      Vprev = Teuchos::null;

      // Get pointer to complete basis
      curind.resize(_curDim);
      for (int i=0; i<_curDim; i++) curind[i] = i;
      Teuchos::RefCountPtr<const MV> curV = MVT::CloneView(*_V,curind);

      // Perform spectral decomposition
      {
        Teuchos::TimeMonitor lcltimer(*_timerCompSF);
        int nevlocal = _curDim;
        int info = _MSUtils.directSolver(_curDim,*_H,0,&S,&_ritzvalues,&nevlocal,10);
        TEST_FOR_EXCEPTION(info != 0,std::logic_error,"Anasazi::BlockKrylovSchur::iterate(): direct solve returned error code.");
        // we did not ask directSolver to perform deflation, so nevLocal 
        TEST_FOR_EXCEPTION(nevlocal != _curDim,std::logic_error,"Anasazi::BlockKrylovSchur::iterate(): direct solve did not compute all eigenvectors."); // this should never happen
      }

      // Sort ritz pairs
      { 
        Teuchos::TimeMonitor lcltimer( *_timerSortEval );

        std::vector<int> _order(_curDim);

        std::vector<ScalarType> _ritzvalues_st(_ritzvalues.size());
        std::copy(_ritzvalues.begin(),_ritzvalues.begin()+_curDim,_ritzvalues_st.begin());

        _sm->sort( this, _curDim, &(_ritzvalues_st[0]), &_order );   // don't catch exception
        
        // Reorder _ritzvalues according to sorting results from _ritzvalues_st
        std::vector<MagnitudeType> _ritzvalues_copy(_ritzvalues);
        for (int i=0; i<_curDim; i++) {
          _ritzvalues[i] = SCT::real(_ritzvalues_copy[_order[i]]);
        }

        // Sort the primitive ritz vectors
        // We need the first _blockSize vectors ordered to generate the next
        // columns immediately below, as well as later, when/if we restart.
        Teuchos::SerialDenseMatrix<int,ScalarType> copyS( S );
        for (int i=0; i<_curDim; i++) {
          blas.COPY(_curDim, copyS[_order[i]], 1, S[i], 1);
        }
      }

      // Create a view matrix of the first _blockSize vectors
      Teuchos::SerialDenseMatrix<int,ScalarType> S1( Teuchos::View, S, _curDim, _blockSize );

      // Compute the new Ritz vectors
      {
	//        Teuchos::TimeMonitor lcltimer( *_timerLocal );
        MVT::MvTimesMatAddMv(ONE,*curV,S1,ZERO,*_X);
      }

      // Update the residual norms
      //_orthman->norm(*_R,&_Rnorms);

      // Update the residual 2-norms 
      //MVT::MvNorm(*_R,&_R2norms);

      // When required, monitor some orthogonalities
      if (_om->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkV = true;
        chk.checkX = true;
        _om->print( Debug, accuracyCheck(chk, ": after local update") );
      }
      else if (_om->isVerbosity( OrthoDetails )) {
        CheckList chk;
        chk.checkX = true;
        _om->print( OrthoDetails, accuracyCheck(chk, ": after local update") );
      }

      // Print information on current iteration
      if (_om->isVerbosity(Debug)) {
        currentStatus( _om->stream(Debug) );
      }
      else if (_om->isVerbosity(IterationDetails)) {
        currentStatus( _om->stream(IterationDetails) );
      }

    } // end while (statusTest == false)

  } // end of iterate()


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Check accuracy, orthogonality, and other debugging stuff
  // 
  // bools specify which tests we want to run (instead of running more than we actually care about)
  //
  // we don't bother checking the following because they are computed explicitly:
  //    H == Prec*R
  //   KH == K*H
  //
  // 
  // checkV : V orthonormal
  //          orthogonal to auxvecs
  // checkX : X orthonormal
  //          orthogonal to auxvecs
  // checkMX: check MX == M*X
  // checkKX: check KX == K*X
  // checkH : H orthonormal 
  //          orthogonal to V and H and auxvecs
  // checkMH: check MH == M*H
  // checkR : check R orthogonal to X
  // checkAux: check that auxiliary vectors are actually orthonormal
  //
  // TODO: 
  //  add checkRitzvalues 
  //
  template <class ScalarType, class MV, class OP>
  std::string BlockKrylovSchur<ScalarType,MV,OP>::accuracyCheck( const CheckList &chk, const string &where ) const 
  {
    stringstream os;
    os.precision(2);
    os.setf(ios::scientific, ios::floatfield);
    MagnitudeType tmp;

    os << " Debugging checks: iteration " << _iter << where << endl;

    // V and friends
    std::vector<int> lclind(_curDim);
    for (int i=0; i<_curDim; i++) lclind[i] = i;
    Teuchos::RefCountPtr<MV> lclV,lclAV;
    lclV = MVT::CloneView(*_V,lclind);
    if (chk.checkV) {
      tmp = _orthman->orthonormError(*lclV);
      os << " >> Error in V^H M V == I  : " << tmp << endl;
      for (unsigned int i=0; i<_auxVecs.size(); i++) {
        tmp = _orthman->orthogError(*lclV,*_auxVecs[i]);
        os << " >> Error in V^H M Q[" << i << "] == 0 : " << tmp << endl;
      }
      /* finish: fix this
      Teuchos::SerialDenseMatrix<int,ScalarType> curH(_curDim,_curDim);
      Teuchos::RefCountPtr<MV> lclAV = MVT::Clone(*_V,_curDim);
      OPT::Apply(*_Op,*lclV,*lclAV);
      MVT::MvTransMv(ONE,*lclV,*lclAV,curH);
      Teuchos::SerialDenseMatrix<int,ScalarType> subH(Teuchos::View,*_H,_curDim,_curDim);
      curH -= subH;
      for (int j=0; j<_curDim; j++) {
        for (int i=j+1; i<_curDim; i++) {
          curH(i,j) = curH(j,i);
        }
      }
      os << " >> Error in V^H A V == H : " << curH.normFrobenius() << endl;
      */
    }

    // X and friends
    // finish: will want some extra tests on X for hermitian problems
    /*
    if (chk.checkX) {
      tmp = _orthman->orthonormError(*_X);
      os << " >> Error in X^H M X == I  : " << tmp << endl;
      for (unsigned int i=0; i<_auxVecs.size(); i++) {
        tmp = _orthman->orthogError(*_X,*_auxVecs[i]);
        os << " >> Error in X^H M Q[" << i << "] == 0 : " << tmp << endl;
      }
    }
    if (chk.checkMX && _hasM) {
      tmp = _MSUtils.errorEquality(_X.get(), _MX.get(), _MOp.get());
      os << " >> Error in MX == M*X     : " << tmp << endl;
    }
    if (chk.checkKX) {
      tmp = _MSUtils.errorEquality(_X.get(), _KX.get(), _Op.get());
      os << " >> Error in KX == K*X     : " << tmp << endl;
    }
    */

    // R: this is not M-orthogonality, but standard euclidean orthogonality
    // finish
    /*
    if (chk.checkR) {
      Teuchos::SerialDenseMatrix<int,ScalarType> xTx(_blockSize,_blockSize);
      MVT::MvTransMv(ONE,*_X,*_R,xTx);
      tmp = xTx.normFrobenius();
      os << " >> Error in X^H R == 0    : " << tmp << endl;
    }
    */

    // Q
    if (chk.checkAux) {
      for (unsigned int i=0; i<_auxVecs.size(); i++) {
        tmp = _orthman->orthonormError(*_auxVecs[i]);
        os << " >> Error in Aux[" << i << "]^H M Aux[" << i << "] == I : " << tmp << endl;
        for (unsigned int j=i+1; j<_auxVecs.size(); j++) {
          tmp = _orthman->orthogError(*_auxVecs[i],*_auxVecs[j]);
          os << " >> Error in Aux[" << i << "]^H M Aux[" << j << "] == 0 : " << tmp << endl;
        }
      }
    }

    os << endl;

    return os.str();
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Print the current status of the solver
  template <class ScalarType, class MV, class OP>
  void 
  BlockKrylovSchur<ScalarType,MV,OP>::currentStatus(ostream &os) 
  {
    os.setf(ios::scientific, ios::floatfield);
    os.precision(6);
    os <<endl;
    os <<"================================================================================" << endl;
    os << endl;
    os <<"                         BlockKrylovSchur Solver Status" << endl;
    os << endl;
    os <<"The solver is "<<(_initialized ? "initialized." : "not initialized.") << endl;
    os <<"The number of iterations performed is " <<_iter<<endl;
    os <<"The block size is         " << _blockSize<<endl;
    os <<"The number of blocks is   " << _numBlocks<<endl;
    os <<"The current basis size is " << _curDim<<endl;
    os <<"The error for the partial Schur decomposition is "<< _schurError <<endl;
    os <<"The number of auxiliary vectors is    " << _numAuxVecs << endl;
    os <<"The number of operations Op*x   is "<<_count_ApplyOp<<endl;

    os.setf(ios_base::right, ios_base::adjustfield);

    /* 
     * FINISH: this will probably have a if (_problem->isHermitian()) thing going on
     * to determine whether to print "Eigenvalue" or "Real/Imag Part"
     *
    if (_initialized) {
      os << endl;
      os <<"CURRENT EIGENVALUE ESTIMATES             "<<endl;
      os << std::setw(20) << "Eigenvalue" 
         << std::setw(20) << "Residual(M)"
         << std::setw(20) << "Residual(2)"
         << endl;
      os <<"--------------------------------------------------------------------------------"<<endl;
      for (int i=0; i<_blockSize; i++) {
        os << std::setw(20) << _ritzvalues[i] 
           << std::setw(20) << _Rnorms[i] 
           << std::setw(20) << _R2norms[i] 
           << endl;
      }
    }
    */

    os <<"================================================================================" << endl;
    os << endl;
  }
  
} // End of namespace Anasazi

#endif

// End of file AnasaziBlockKrylovSchur.hpp
