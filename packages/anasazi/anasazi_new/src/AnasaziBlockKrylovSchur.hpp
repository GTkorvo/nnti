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

#ifndef ANASAZI_BLOCK_ARNOLDI_HPP
#define ANASAZI_BLOCK_ARNOLDI_HPP

#include "AnasaziEigensolver.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziSortManager.hpp"
#include "AnasaziOutputManager.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"

/*!	\class Anasazi::BlockKrylovSchur

	\brief This class implements the Restarted Block Krylov Schur Method,
	an iterative method for solving eigenvalue problems.

	\author Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {
  
  template <class TYPE>
  class BlockKrylovSchur : public Eigensolver<TYPE> { 
  public:
    //@{ \name Constructor/Destructor.
    
    //! %Anasazi::BlockKrylovSchur constructor.
    BlockKrylovSchur( Eigenproblem<TYPE>& problem, 
		  SortManager<TYPE>& sm,
		  OutputManager<TYPE>& om,
		  const TYPE tol=1.0e-6,
		  const int length=25, 
		  const int step=25, 
		  const int restarts=0 
		  );
    
    //! %Anasazi::BlockKrylovSchur destructor.
    virtual ~BlockKrylovSchur();
    //@}
    
    //@{ \name Solver application methods.
    
    /*! \brief This method performs a given number of steps of the BlockArnoldi
      Method, returning upon completion or convergence.
    */
    void iterate( const int steps=1 );
    
    /*! \brief This method uses iterate to compute approximate solutions to the
      original problem.  It may return without converging if it has taken the
      maximum number of iterations or numerical breakdown is observed.
    */
    void solve();
    //@}
    
    //@{ \name Solution return methods.
    
    //! This method returns the computed Ritz values.
    const TYPE * GetRitzValues() { return(_ritzvalues); };
    
    //! This method returns the Ritz residuals for the computed eigenpairs.
    const TYPE * GetRitzResiduals() { return(_ritzresiduals); };

    //! Get the current iteration count.
    int GetNumIters() const { return(_iter); };
    
    //! Get the current restart count of the iteration method.
    /*! Some eigensolvers can perform restarts (i.e.Arnoldi) to reduce memory
      and orthogonalization costs.  For other eigensolvers that don't
      perform restarts (i.e. LOBPCG), this is not a valid stopping criteria.
    */
    int GetNumRestarts() const { return(_restarts); };
    
    //! Get the solvers native residuals for the current eigenpairs. 
    /*! This is not be the same as the true residuals for most solvers. Sometimes the native
      residuals are not in multivector form, so the norm type is solver dependent.  
      
      \note
      <ol>
      <li> If the native residual is in multivector form then a non-null pointer will be
      returned, else the normvec will be populated with the current residual norms. 
      <li> If the native residual is returned in multivector form, the memory is managed
      by the calling routine.
      </ol>
    */
    MultiVec<TYPE>* GetNativeResiduals( TYPE* normvec ) const {};
    
    /*! \brief Get a constant reference to the current linear problem, 
      which may include a current solution.
    */
    Eigenproblem<TYPE>& GetEigenproblem() const { return(_problem); };
    
    //@}
    
    //@{ \name Output methods.
    
    //! This method requests that the solver print out its current status to screen.
    void currentStatus();
    //@}
  private:
    void QRFactorization( MultiVec<TYPE>&, Teuchos::SerialDenseMatrix<int,TYPE>& );
    void ComputeSchurForm( const bool apply );
    void SortSchurForm( Teuchos::SerialDenseMatrix<int,TYPE>& H, Teuchos::SerialDenseMatrix<int,TYPE>& Q );
    void BlockReduction();
    void BlkOrth( MultiVec<TYPE>& Vec_in, const int j );
    void BlkOrthSing( MultiVec<TYPE>& Vec_in, const int j );
    void ComputeEvecs();
    void Restart();
    void SetBlkTols();
    void CheckBlkArnRed( const int j );
    void CheckSchurVecs( const int j ); 

    Eigenproblem<TYPE> &_problem; // must be passed in by the user
    SortManager<TYPE> &_sm; // must be passed in by the user
    OutputManager<TYPE> &_om; // must be passed in by the user

    MultiVec<TYPE> *_basisvecs;
    Teuchos::SerialDenseMatrix<int,TYPE> _hessmatrix;
    const int _length, _restarts, _step;
    const TYPE _residual_tolerance;
    TYPE *_ritzresiduals, *_ritzvalues, *_ritzvaluesi;
    int *_order;
    int _restartiter, _iter, _jstart, _jend, _nevblock;
    int _offset, _maxoffset;
    bool _isdecompcurrent, _isevecscurrent, _exit_flg, _dep_flg;
    TYPE _schurerror, _dep_tol, _blk_tol, _sing_tol, _def_tol;

    // Information obtained from the eigenproblem
    Operator<TYPE> *_Op;
    Operator<TYPE> *_B;
    MultiVec<TYPE> *_evecs;
    const int _nev, _block;  
    TYPE *_evals;
  };
  //
  // Implementation
  //
  // Note: I should define a copy constructor and overload = because of the use of new
  //
  template <class TYPE>
  BlockKrylovSchur<TYPE>::BlockKrylovSchur(Eigenproblem<TYPE> & problem, 
				   SortManager<TYPE> & sm,
				   OutputManager<TYPE> & om,
				   const TYPE tol, 
				   const int length, 
				   const int step, 
				   const int restarts
				   ): 
    _problem(problem), 
    _sm(sm),
    _om(om),
    _basisvecs(0), 
    _hessmatrix(),
    _length(length), 
    _restarts(restarts),
    _step(step),
    _residual_tolerance(tol),
    _ritzresiduals(0), 
    _ritzvalues(0), 
    _ritzvaluesi(0),
    _order(0),
    _restartiter(0), 
    _iter(0), 
    _jstart(0), 
    _jend(0), 
    _nevblock(0),
    _offset(0),
    _maxoffset(0),
    _isdecompcurrent(false),
    _isevecscurrent(false),
    _exit_flg(false),
    _dep_flg(false),
    _schurerror(1.0), 
    _dep_tol(1.0), 
    _blk_tol(1.0),
    _sing_tol(1.0),
    _def_tol(1.0),
    _Op(problem.GetOperator()),
    _B(problem.GetB()),
    _evecs(problem.GetEvecs()), 
    _nev(problem.GetNEV()), 
    _block(problem.GetBlockSize()), 
    _evals(problem.GetEvals()) 
    {     
    //
    // Determine _nevblock : how many blocks it will take to contain the _nev eigenvalues/vectors
    // NOTE: An additional block is kept if _nev is a perfect multiple of _block because of the
    // potential presence of complex eigenvalue pairs.  Additional blocks can be retained, up to
    // _maxoffset if the block ends with one eigenvalue of a complex conjugate pair.
    //
    _nevblock = _nev/_block + 1;
    // TEST CODE:  This part has changed from saving another block to compare with ARPACK.
    //_nevblock = _nev/_block;  
    //if (_nev%_block) 
    //_nevblock++;    
    _maxoffset = (_length-_nevblock)/2;
    //
    // Retrieve the initial vector from the Anasazi::Eigenproblem.
    //
    MultiVec<TYPE>* ivec = _problem.GetInitVec();
    assert(ivec!=NULL);
    //
    assert(_length>0); assert(_step>0);
    //
    // Make room for theArnoldi vectors and F.
    //
    _basisvecs = ivec->Clone((_length+1)*_block); assert(_basisvecs!=NULL);
    //
    // Create the rectangular Hessenberg matrix
    //
    _hessmatrix.shape((_length+1)*_block, _length*_block); 
    //
    // Create the vectors for eigenvalues and their residual errors and
    // initialize them.
    //
    _ritzvalues = new TYPE[ 2* (_block*_length) ]; assert(_ritzvalues!=NULL);  
    _ritzvaluesi = _ritzvalues + _block*_length;
    _ritzresiduals = new TYPE[ _block*_length ]; assert(_ritzresiduals!=NULL);
    _order = new int[ _block*_length ]; assert(_order!=NULL);
    const TYPE one = 1.0, zero = 0.0;
    for (int i=0; i< _block*_length; i++) {
      _ritzvalues[i] = zero; _ritzvaluesi[i] = zero;
      _ritzresiduals[i] = one;
    }			
    //
    //  Set the tolerances for block orthogonality
    //
    SetBlkTols();  
  }
  
  template <class TYPE>
  void BlockKrylovSchur<TYPE>::SetBlkTols() {
    const TYPE two = 2.0;
    TYPE eps;
    char precision = 'P';
    Teuchos::LAPACK<int,TYPE> lapack;
    eps = lapack.LAMCH(precision);
    _blk_tol = 10*sqrt(eps);
    _sing_tol = 10 * eps;
    _dep_tol = 1/sqrt(two);
    _def_tol = eps;
  }
  
  template <class TYPE>
  BlockKrylovSchur<TYPE>::~BlockKrylovSchur() 
  {
    if (_basisvecs) delete _basisvecs;
    if (_ritzresiduals) delete [] _ritzresiduals;
    if (_ritzvalues) delete [] _ritzvalues;
    if (_order) delete [] _order;
  }
  
  template <class TYPE>
  void BlockKrylovSchur<TYPE>::currentStatus() {
    int i;
    if (_om.doOutput(-1)) {
      cout<<" "<<endl;
      cout<<"********************CURRENT STATUS********************"<<endl;
      cout<<"Iterations :\t"<<_iter<<endl;
      
      if (_restartiter > _restarts) 
	cout<<"Restarts :\t"<<_restartiter-1<<" of\t"<< _restarts<<endl;
      else
	cout<<"Restarts :\t"<<_restartiter<<" of\t"<< _restarts<<endl;
      
      cout<<"Block Size :\t"<<_block<<endl;
      cout<<"Requested Eigenvalues : "<<_nev<<endl;
      cout<<"Residual Tolerance : "<<_residual_tolerance<<endl;	
      cout<<"Error for the partial Schur decomposition is : "<< _schurerror <<endl;
      //
      //  Determine status of solver and output information correctly.
      //
      if ( _schurerror < _residual_tolerance ) {
	cout<<"------------------------------------------------------"<<endl;
	cout<<"Computed Eigenvalues: "<<endl;
      } else {
	if (_exit_flg && _iter != _length+_restarts*(_length-_nevblock)) {
	  cout<<"ERROR: Complete orthogonal basis could not be computed"<<endl;
	}
	cout<<"------------------------------------------------------"<<endl;
	cout<<"Current Eigenvalue Estimates: "<<endl;
      }
      //
      //  Print out current computed eigenvalues.  If we don't have all the requested
      //  eigenvalues yet, print out the ones we have.
      //
      int _nevtemp = _nev;
      if (_jstart < _nevblock) { _nevtemp = _jstart*_block; }
      //
      if (_problem.IsSymmetric()) {
	cout<<"Eigenvalue\tRitz Residual"<<endl;
	cout<<"------------------------------------------------------"<<endl;
	if ( _nevtemp == 0 ) {
	  cout<<"[none computed]"<<endl;
	} else {
	  for (i=0; i<_nevtemp; i++) {
	    cout.width(10);
	    cout<<_evals[i]<<"\t"<<_ritzresiduals[i]<<endl;
	  }
	}
	cout<<"------------------------------------------------------"<<endl;
      } else {
	cout<<"Real Part\tImag Part\tRitz Residual"<<endl;
	cout<<"------------------------------------------------------"<<endl;
	if ( _nevtemp == 0 ) {
	  cout<<"[none computed]"<<endl;
	} else {
	  for (i=0; i<_nevtemp; i++) {
	    cout.width(10);
	    cout<<_evals[i]<<"\t"<<_evals[_nev+i]<<"\t\t"<<_ritzresiduals[i]<<endl;
	  }
	}
	cout<<"------------------------------------------------------"<<endl;
	cout<<" "<<endl;
      }
      cout<<"******************************************************"<<endl;
    }	
  }
  
  template <class TYPE>
  void BlockKrylovSchur<TYPE>::solve () {
    //int rem_iters = _length+_restarts*(_length-_nevblock)-_iter;
    //
    // Right now the solver will just go the remaining iterations, but this design will allow
    // for checking of the residuals every so many iterations, independent of restarts.
    //
    while (_schurerror > _residual_tolerance && _restartiter <= _restarts && !_exit_flg) {
      iterate( _step );
    }
    //
    // Compute the current approximate eigenvectors before returning.
    //
    ComputeEvecs();    
    //
  }
  
  template <class TYPE>
  void BlockKrylovSchur<TYPE>::iterate(const int steps) {
    int i,temp;
    int tempsteps = steps;
    //
    // If this is the first steps of Block Krylov Schur, initialize the first block of _basisvecs
    //
    if (!_iter) {
      int *index = new int[ _block ]; assert(index!=NULL);
      for (i=0; i<_block; i++) {
	index[i] = i;
      }
      //
      // Copy the first _block of the initial vectors into the first _block
      // of _basisvecs, any additional vectors will be ignored.
      //
      _basisvecs->SetBlock( *(_problem.GetInitVec()), index, _block );
      //
      // Orthogonalize the first block of vectors.
      //      
      MultiVec<TYPE>* U_vec = _basisvecs->CloneView( index,_block );
      assert(U_vec!=NULL);
      Teuchos::SerialDenseMatrix<int,TYPE> G10( _block,_block );
      QRFactorization( *U_vec, G10 );
      delete U_vec;
      delete [] index;
    }				
    //
    // Leave the iteration method now if the orthogonal subspace can't be extended.
    //
    if (_exit_flg) { return; }	
    //			
    // Now we go the number of steps requested by the user.  This may cause
    // a restart or hit the number of maximum iterations (restarts).  
    //
    while(tempsteps > 0 && _restartiter <= _restarts && !_exit_flg) {
      _isevecscurrent = false;
      // If we don't need to restart, just get it over with and return.
      if (_jstart+tempsteps < _length) {
	_jend = _jstart+tempsteps;
	BlockReduction();
	//
	// We need to leave before we move the pointer if the orthogonalization failed.
	if (_exit_flg ) { break; } 
	//
	// Move the pointer and update the iteration count.
	//
	_iter += tempsteps;
	tempsteps = 0;
	_jstart = _jend;  
	ComputeSchurForm( false );		
	_isdecompcurrent = false;
	// Output current information if necessary
	if (_om.doOutput(0)) {
	  currentStatus();
	}
      }
      // Finish off this factorization and restart.
      else {  
	_jend = _length;
	BlockReduction();
	//
	// We need to leave before we move the pointer if the orthogonalization failed.
	if (_exit_flg ) { break; } 
	//
	// Move the pointer and update the iteration count.
	//
	temp = _length-_jstart;
	_iter += temp;
	tempsteps -= temp;
	_jstart = _length; 
	//
	//  Compute the Schur factorization and prepare for a restart.  Don't
	//  compute restart if at end of iterations.
	//
	if (_restartiter < _restarts) {
	  ComputeSchurForm( true );  
	  Restart();  
	  _isdecompcurrent = true;
	  _restartiter++;
	} else {
	  ComputeSchurForm( false );
	  _restartiter++;
	  _isdecompcurrent = false;
	}
	// Output current information if necessary
	if (_om.doOutput(0)) {
	  currentStatus();
	}
      }
    }
  }
  
  
  template<class TYPE>
  void BlockKrylovSchur<TYPE>::BlockReduction () {
    int i,j;
    ReturnType ret;
    int *index = new int[ _block ]; assert(index!=NULL);
    MultiVec<TYPE> *U_vec=0, *F_vec=0;
    
    for ( j = _jstart; j < _jend; j++ ) {
      //
      // Associate the j-th block of _basisvecs with U_vec.
      //
      for ( i=0; i<_block; i++ ) {
	index[i] = j*_block+i;
      }
      U_vec = _basisvecs->CloneView(index, _block);
      assert(U_vec!=NULL);
      //
      // Associate (j+1)-st block of ArnoldiVecs with F_vec.
      //
      //for ( i=0; i<_block; i++ ) {
      //	index[i] = (j+1)*_block+i;
      //}
      F_vec = _basisvecs->Clone(_block);
      //F_vec = _basisvecs->CloneView(index, _block);
      //assert(F_vec!=NULL);
      //
      //  Compute F_vec = OP * U_vec
      //
      ret =_Op->Apply( *U_vec, *F_vec ); 
      //
      // Use previous dependency information to decide which orthogonalization
      // method to use for the new block.  The global variable _dep_flg tells us
      // if we've had problems with orthogonality before.  If no problems have
      // been detected before we will use standard block orthogonalization.
      // However, if there are problems with this, we will use a more stringent
      // approach.
      //
      if (!_dep_flg) {
	BlkOrth( *F_vec, j );
      }
      //
      // If any block dependency was detected previously, then the more stringent
      // orthogonalization will be used.  If this method can't resolve the
      // dependency, then the _exit_flg will be set indicating that we can't proceed
      // any further.
      //			
      if (_dep_flg) {
	BlkOrthSing( *F_vec, j );
      }
      //
      delete U_vec; U_vec=0;
      delete F_vec; F_vec=0;
      //
      // If we cannot go any further with the factorization, then we need to exit
      // this method.
      //
      if (_exit_flg) { break; }
    }
    delete [] index;
  } // end BlockReduction()
  
  
  template<class TYPE>
  void BlockKrylovSchur<TYPE>::BlkOrth( MultiVec<TYPE>& Vec_in, const int j ) {
    //
    // Orthogonalization is first done between the new block of
    // vectors and all previous blocks, then the vectors within the
    // new block are orthogonalized.
    //
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    const int max_num_orth = 2;
    int i, k, row_offset, col_offset;
    int * index = new int[ (_length+1)*_block ]; assert(index!=NULL);
    TYPE * norm1 = new TYPE[_block]; assert(norm1!=NULL);
    TYPE * norm2 = new TYPE[_block]; assert(norm2!=NULL);
    ReturnType ret; 
    //
    // Associate (j+1)-st block of ArnoldiVecs with F_vec.
    //
    for ( i=0; i<_block; i++ ) {
      index[i] = (j+1)*_block+i;
    }
    MultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _block);
    assert(F_vec!=NULL);
    F_vec->MvAddMv(one, Vec_in, zero, Vec_in);
    //
    // Zero out the full block column of the Hessenberg matrix
    // even though we're only going to set the coefficients in
    // rows [0:(j+1)*_block-1]
    //
    int n_row = _hessmatrix.numRows();
    //
    for ( k=0; k<_block; k++ ) {
      for ( i=0; i<n_row ; i++ ) {
	_hessmatrix( i, j*_block+k ) = zero;
      }
    }
    //
    // Grab all previous Arnoldi vectors
    //
    int num_prev = (j+1)*_block;
    for (i=0; i<num_prev; i++){
      index[i] = i;
    }
    MultiVec<TYPE>* V_prev = _basisvecs->CloneView(index,num_prev);
    assert(V_prev!=NULL);
    //
    // Create a matrix to store the product trans(V_prev)*B*F_vec
    //
    Teuchos::SerialDenseMatrix<int,TYPE> dense_mat(num_prev, _block );
    //
    F_vec->MvNorm(norm1);
    //
    // Check the norm of the candidate block of vectors to make sure they're
    // not zero.  [ This might happen if the matrix is the zero matrix ]
    //
    for (i=0; i<_block; i++) {
      if (norm1[i] == zero) {
	_dep_flg = true;
	if (_om.doOutput(2)){
	  cout << "Col " << num_prev+i << " is the zero vector" << endl;
	  cout << endl;
	}
      }	  
    }
    //
    // Perform two steps of block classical Gram-Schmidt so that
    // F_vec is B-orthogonal to the columns of V_prev.
    //
    for ( int num_orth=0; num_orth<max_num_orth; num_orth++ ) {
      //
      // Compute trans(V_prev)*B*F_vec and store in the j'th diagonal
      // block of the Hessenberg matrix
      //
      ret = _problem.InnerProd( *V_prev, *F_vec, dense_mat );
      //
      // Update the orthogonalization coefficients for the j-th block
      // column of the Hessenberg matrix.
      //
      for ( k=0; k<_block; k++ ) {
	for ( i=0; i<num_prev; i++ ) {
	  _hessmatrix( i, j*_block+k ) += dense_mat(i,k);
	}
      }
      //
      // F_vec <- F_vec - V(0:(j+1)*block-1,:) * H(0:num_prev-1,j:num_prev-1)
      //
      F_vec->MvTimesMatAddMv( -one, *V_prev, dense_mat, one );
      //
    } // end for num_orth=0;...)
    //
    F_vec->MvNorm(norm2);
    //
    // Check to make sure the new block of Arnoldi vectors are
    // not dependent on previous Arnoldi vectors.  
    //
    for (i=0; i<_block; i++){
      if (norm2[i] < norm1[i] * _blk_tol) {
	_dep_flg = true;
	if (_om.doOutput(2)){
	  cout << "Col " << num_prev+i << " is dependent on previous "
	       << "Arnoldi vectors in V_prev" << endl;
	  cout << endl;
	}
      }
    } // end for (i=0;...)
    //
    if (_om.doOutput(2)) {
      CheckBlkArnRed(j);
    }
    //
    // If dependencies have not already been detected, compute
    // the QR factorization of the next block. Otherwise,
    // this block of Arnoldi vectors will be re-computed via and
    // implementation of A. Ruhe's block Arnoldi.
    //
    if (!_dep_flg) {
      //
      // Compute the QR factorization of F_vec
      //
      row_offset = (j+1)*_block; col_offset = j*_block;
      Teuchos::SerialDenseMatrix<int,TYPE> sub_block_hess(Teuchos::View, _hessmatrix, _block, _block, 
							  row_offset, col_offset);
      
      QRFactorization( *F_vec, sub_block_hess );
    }
    //
    delete F_vec;
    delete V_prev;
    delete [] index;
    delete [] norm1;
    delete [] norm2;
    //
  }  // end BlkOrth()
  
  
  template<class TYPE>
  void BlockKrylovSchur<TYPE>::BlkOrthSing( MultiVec<TYPE>& Vec_in, const int j ) {
    //
    // This is a variant of A. Ruhe's block Arnoldi
    // The orthogonalization of the vectors F_vec is done
    // one at a time. If a dependency is detected, a random
    // vector is added and orthogonalized against all previous
    // Arnoldi vectors.
    //
    const int IntOne = 1;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    int i, k, num_prev;
    int * index = new int[ (_length+1)*_block ]; assert(index!=NULL);
    Teuchos::SerialDenseVector<int,TYPE> dense_vec;
    TYPE norm1[IntOne];
    TYPE norm2[IntOne];
    ReturnType ret;
    //
    // Place the candidate vectors Vec_in into the (j+1)-st block of ArnoldiVecs.
    //
    for ( i=0; i<_block; i++ ) {
      index[i] = (j+1)*_block+i;
    }
    _basisvecs->SetBlock( Vec_in, index, _block ); 
    //
    // Zero out the full block column of the Hessenberg matrix
    //
    int n_row = _hessmatrix.numRows();
    //
    for ( k=0; k<_block; k++ ) {
      for ( i=0; i<n_row ; i++ ) {
	_hessmatrix(i, j*_block+k) = zero;
      }
    }
    //
    MultiVec<TYPE> *q_vec=0, *Q_vec=0, *tptr=0;
    tptr = _basisvecs->Clone(IntOne); assert(tptr!=NULL);
    //
    // Start a loop to orthogonalize each of the _block
    // columns of the (j+1)-st block of _basisvecs against all 
    // the others.
    //
    for (int iter=0; iter<_block; iter++){
      num_prev = (j+1)*_block + iter; // number of previous _basisvecs
      dense_vec.size(num_prev);
      //
      // Grab the next column of _basisvecs
      //
      index[0] = num_prev;
      q_vec = _basisvecs->CloneView(index, IntOne); assert(q_vec!=NULL);
      //
      // Grab all previous columns of _basisvecs
      //
      for (i=0; i<num_prev; i++){
	index[i] = i;
      }
      Q_vec = _basisvecs->CloneView(index, num_prev); assert(Q_vec!=NULL);
      //
      // Create matrix to store product trans(Q_vec)*B*q_vec
      //
      // Do one step of classical Gram-Schmidt B-orthogonalization
      // with a 2nd correction step if needed.
      //
      q_vec->MvNorm(norm1);
      //
      // Leave if this is the zero vector, there is no more we can do here.
      //
      if (norm1[0] == zero) { 
	if (_om.doOutput(2)) {
	  cout << "Column " << num_prev << " of _basisvecs is the zero vector" 
	       << endl<<endl;
	}
	_exit_flg = true; 
	break; 
      }
      //
      // Compute trans(Q_vec)*B*q_vec
      //
      ret = _problem.InnerProd( *Q_vec, *q_vec, dense_vec );
      //
      // Sum results [0:num_prev-1] into column (num_prev-_block)
      // of the Hessenberg matrix
      //
      for (k=0; k<num_prev; k++){
	_hessmatrix(k, j*_block+iter) += dense_vec(k);
      }
      //
      // Compute q_vec<- q_vec - Q_vec * dense_vec
      //
      q_vec->MvTimesMatAddMv(-one, *Q_vec, dense_vec, one);
      //
      q_vec->MvNorm(norm2);
      //
      if (norm2[0] < norm1[0] * _dep_tol) {
	//
	// Repeat process with newly computed q_vec
	//
	// Compute trans(Q_vec)*q_vec
	//
	ret = _problem.InnerProd( *Q_vec, *q_vec, dense_vec );
	//
	// Sum results [0:num_prev-1] into column (num_prev-_block)
	// of the Hessenberg matrix
	//
	for (k=0; k<num_prev; k++){
	  _hessmatrix(k, j*_block+iter) += dense_vec(k);
	}
	//
	// Compute q_vec<- q_vec - Q_vec * dense_vec
	//
	q_vec->MvTimesMatAddMv(-one, *Q_vec, dense_vec, one);
	//
	q_vec->MvNorm(norm2);
      }
      //
      // Check for linear dependence
      //
      if (norm2[0] < norm1[0] * _sing_tol) {
	if (_om.doOutput(2)) {
	  cout << "Column " << num_prev << " of _basisvecs is dependent" 
	       << endl<<endl;
	}
	//
	// Create a random vector and orthogonalize it against all
	// previous cols of _basisvecs
	// We could try adding a random unit vector instead -- not
	// sure if this would make any difference.
	//
	tptr->MvRandom();
	tptr->MvNorm(norm1);
	//
	// This code  is automatically doing 2 steps of B-orthogonalization
	// after adding a random vector. We could do one step of
	// orthogonalization with a correction step if needed.
	//
	for (int num_orth=0; num_orth<2; num_orth++){
	  ret = _problem.InnerProd( *Q_vec, *tptr, dense_vec );
	  // Note that we don't change the entries of the
	  // Hessenberg matrix when we orthogonalize a
	  // random vector
	  tptr->MvTimesMatAddMv(-one, *Q_vec, dense_vec, one);
	}
	//
	tptr->MvNorm(norm2);
	//
	if (norm2[0] > norm1[0] * _sing_tol){
	  // Copy vector into the current column of _basisvecs
	  q_vec->MvAddMv( one, *tptr, zero, *tptr );
	  q_vec->MvNorm(norm2);
	  //
	  // Normalize the new q_vec
	  //
	  TYPE rjj = one/norm2[0];
	  q_vec->MvAddMv( rjj, *q_vec, zero, *q_vec );
	  //
	  // Enter a zero in the [(j+1)*_block + iter] row in the
	  // [(j*_block + iter] column of the Hessenberg matrix
	  //
	  _hessmatrix((j+1)*_block+iter, j*_block+iter) = zero;
	}
	else {
	  // Can't produce a new orthonormal basis vector
	  // Clean up and exit this block Arnoldi factorization!
	  _exit_flg = true;
	  delete [] index;
	  delete q_vec; q_vec=0;
	  delete Q_vec; Q_vec=0;
	  delete tptr; tptr=0;
	  return;
	}
      }
      else {
	//
	// Normalize the new q_vec
	//
	TYPE rjj = one/norm2[0];
	q_vec->MvAddMv( rjj, *q_vec, zero, *q_vec );
	//
	// Enter norm of q_vec to the [(j+1)*_block + iter] row
	// in the [(j*_block + iter] column of the Hessenberg matrix
	//
	_hessmatrix((j+1)*_block+iter, j*_block+iter) = norm2[0];
      } // end else ...
    } // end for (i=0;...)
    //
    if (_om.doOutput(2)){
      cout << "Checking Orthogonality after BlkOrthSing()"
	   << " Iteration: " << j << endl<<endl;
      CheckBlkArnRed(j);
    }
    //
    //      free heap space
    //
    delete [] index;
    delete q_vec; q_vec=0;
    delete Q_vec; Q_vec=0;
    delete tptr; tptr=0;
  } // end BlkOrthSing()

  
  template<class TYPE>
  void BlockKrylovSchur<TYPE>::QRFactorization (MultiVec<TYPE>& VecIn, 
					    Teuchos::SerialDenseMatrix<int,TYPE>& FouierR) {
    int i,j,k;
    int nb = VecIn.GetNumberVecs(); assert (nb == _block);
    int *index = new int[nb]; assert(index!=NULL);
    const int IntOne=1;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    bool addvec = false, flg = false;
    ReturnType ret;
    //
    TYPE norm1[IntOne];
    TYPE norm2[IntOne];
    MultiVec<TYPE> *qj = 0, *Qj = 0, *tptr = 0;
    tptr = _basisvecs->Clone(IntOne); assert(tptr!=NULL);
    //
    // Zero out the array that will contain the Fourier coefficients.
    //
    for ( j=0; j<nb; j++ ) {
      for ( i=0; i<nb; i++ ) {
	FouierR(i,j) = zero;
      }
    }
    //
    // Start the loop to orthogonalize the nb columns of VecIn.
    //
    for ( j=0; j<nb; j++ ) {
      //
      flg = false;
      //
      // Grab the j-th column of VecIn (the first column is indexed to 
      // be the zero-th one).
      //
      index[0] = j;
      qj = VecIn.CloneView(index, IntOne); assert(qj!=NULL);
      //
      // If we are beyong the 1st column, orthogonalize against the previous
      // vectors in the current block.
      //
      if ( j ) {
	for ( i=0; i<j; i++ ) {
	  index[i] = i;
	}
	//
	// Grab the first j columns of VecIn (that are now an orthogonal
	// basis for first j columns of the entering VecIn).
	//
	Qj = VecIn.CloneView(index, j);
	Teuchos::SerialDenseVector<int,TYPE> rj(j);
	_problem.MvNorm( *qj, norm1 );
	//
	// Do one step of classical Gram-Schmidt orthogonalization
	// with a second correction step if needed
	//
	// Determine the Fouier coefficients for B-orthogonalizing column
	// j of VecIn against columns 0:j-1 of VecIn. In other words,
	// result = trans(Qj)*B*qj.
	//
	ret = _problem.InnerProd( *Qj, *qj, rj );
	//
	// Sum results[0:j-1] into column j of R.
	//
	for ( k=0; k<j; k++ ) {
	  FouierR(k,j) += rj(k);
	}
	//
	// Compute qj <- qj - Qj * rj.
	//
	qj->MvTimesMatAddMv(-one, *Qj, rj, one);
	//
	_problem.MvNorm( *qj, norm2 );			
	//
	if (norm2[0] < norm1[0] * _dep_tol){
	  //
	  // Repeat process with newly computed qj
	  //
	  ret = _problem.InnerProd( *Qj, *qj, rj );
	  //    				
	  // Sum results[0:j-1] into column j of R.
	  //
	  for ( k=0; k<j; k++ ) {
	    FouierR(k,j) += rj(k);
	  }
	  //
	  // Compute qj <- qj - Qj * rj.
	  //
	  qj->MvTimesMatAddMv(-one, *Qj, rj, one);
	  //
	  _problem.MvNorm( *qj, norm2 );
	}
	//
	// Check for dependencies
	//
	if (_iter) {
	  // This is not the 1st block. A looser tolerance is used to
	  // determine dependencies. If a dependency is detected, a flag
	  // is set so we can back out this method and out of BlkOrth.
	  // The method BlkOrthSing is used to construct the new block
	  // of orthonormal basis vectors one at a time. If a dependency
	  // is detected within this method, a random vector is added
	  // and orthogonalized against all previous basis vectors.
	  //
	  if (norm2[0] < norm1[0] * _blk_tol) {
	    if (_om.doOutput(2)) {
	      cout << "Column " << j << " of current block is dependent"<<endl;
	    }
	    _dep_flg = true;
	    delete qj; delete Qj; delete tptr;
	    delete [] index;
	    return;
	  }
	}
	else {
	  // This is the 1st block of basis vectors.
	  // Use a tighter tolerance to determine dependencies, because
	  // if a dependency is detected we will be adding a random
	  // vector and orthogonalizing it against previous vectors
	  // in the 1st block
	  //
	  if (norm2[0] < norm1[0] * _sing_tol) {
	    // The 1st block of vectors are dependent
	    // Add a random vector and orthogonalize it against
	    // previous vectors in block.
	    //
	    addvec = true;
	    Teuchos::SerialDenseVector<int,TYPE> tj(j);
	    //
	    tptr->MvRandom();
	    _problem.MvNorm( *tptr, norm1 );
	    //
	    for (int num_orth=0; num_orth<2; num_orth++){
	      ret = _problem.InnerProd( *Qj, *tptr, tj );
	      tptr->MvTimesMatAddMv(-one, *Qj, tj, one);
	    }
	    _problem.MvNorm( *tptr, norm2 );
	    //
	    if (norm2[0] > norm1[0] * _sing_tol){
	      // Copy vector into current column of _basisvecs
	      qj->MvAddMv(one, *tptr, zero, *tptr);
	    }
	    else {
	      _exit_flg = true;
	      delete qj; delete Qj; delete tptr;
	      delete [] index;
	      return;
	    }
	  }
	} // if (_iter) ...
      } // if (j) ...
      //
      // If we have not exited, compute the norm of column j of
      // VecIn (qj), then normalize qj to make it into a unit vector
      //
      TYPE normq[IntOne];
      _problem.MvNorm( *qj, normq );
      //
      TYPE rjj = one / normq[0];
      qj->MvAddMv ( rjj, *qj, zero, *qj );
      //
      if (addvec){
	// We've added a random vector, so
	// enter a zero in j'th diagonal element of R
	FouierR(j,j) = zero;
      }
      else {
	FouierR(j,j) = normq[0];
      }
      delete qj; delete Qj;
    } // for (j=0; j<nb; j++) ...
    //
    delete tptr;
    delete [] index;	
  }
  
  template<class TYPE>
  void BlockKrylovSchur<TYPE>::ComputeEvecs() {
    //
    int i=0,j=0,k=0;
    int n=_jstart*_block, info=0;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    Teuchos::LAPACK<int,TYPE> lapack;
    Teuchos::BLAS<int,TYPE> blas;
    MultiVec<TYPE>* basistemp=0;
    Teuchos::SerialDenseMatrix<int,TYPE> Q(n,n);
    int * index = new int [ n ]; assert(index!=NULL);
    //
    // Initialize the eigenvectors to zero.
    //
    _evecs->MvInit( 0.0 );
    //
    //  Set the index array.
    //
    for (k=0; k<n; k++) { index[k] = k; }
    //
    //  By some chance the number of eigenvalues requested may be more than
    //  the size of the factorization.  To prevent a seg fault in creating the
    //  eigenvectors, determine if this is the current situation.
    //
    int curr_nev = _nev;
    if ( n < _nev ) { curr_nev = n; }
    //
    //  Get a view into the current Hessenberg matrix.
    //
    Teuchos::SerialDenseMatrix<int,TYPE> Hj(Teuchos::Copy, _hessmatrix, n, n, i, j);
    //
    //  If the Krylov-Schur decomposition is not current, compute residuals
    //  like we are going to restart to update decomposition.
    //
    if (!_isdecompcurrent) { 
      //
      // Update the Krylov basis using the Schur form.  
      // Take into account that deflated blocks need not be updated.
      //	
      SortSchurForm( Hj, Q );
      basistemp = _basisvecs->Clone( n );
      MultiVec<TYPE>* basistemp2 = _basisvecs->CloneView( index, n );
      basistemp->MvTimesMatAddMv ( one, *basistemp2, Q, zero );
      delete basistemp2;
    } else {
      //
      // We can aquire the Ritz vectors from the current decomposition.
      //
      basistemp = _basisvecs->CloneCopy( index, n );
    }
    //
    // Check the Schur form.
    //
    if (_om.doOutput(2))
      CheckSchurVecs( _jstart );
    //
    //  If the operator is symmetric, then the Ritz vectors are the eigenvectors.
    //  So, copy the Ritz vectors.  Else, we need to compute the eigenvectors of the
    //  Schur form to compute the eigenvectors of the non-symmetric operator.
    //
    if (_problem.IsSymmetric()) {
      _evecs->SetBlock( *basistemp, index, curr_nev );
    } else {  
      //
      //  Now compute the eigenvectors of the Schur form
      //  Reset the dense matrix and compute the eigenvalues of the Schur form.
      //
      int lwork = 4*n;
      TYPE *work = new TYPE[lwork]; assert(work!=NULL);
      int *select = new int[ n ];	  
      char * side = "R";
      char * howmny = "A";
      int mm, ldvl = 1;
      TYPE *vl = new TYPE[ ldvl ];
      lapack.TREVC( *side, *howmny, select, n, Hj.values(), Hj.stride(), vl, ldvl,
		    Q.values(), Q.stride(), n, &mm, work, &info );
      assert(info==0);
      delete [] work;
      delete [] select;
      delete [] vl;
      //
      //  Convert back to approximate eigenvectors of the operator and compute their norms.
      //
      TYPE* evecnrm = new double[ n ];
      MultiVec<TYPE>* evecstemp = _basisvecs->Clone( n );
      evecstemp->MvTimesMatAddMv( one, *basistemp, Q, zero );
      evecstemp->MvNorm( evecnrm );
      //
      // Sort the eigenvectors.
      //
      int conjprs=0;
      int * indexi = new int [ curr_nev+1 ];
      MultiVec<TYPE> *evecstempr, *evecr1;
      TYPE t_evecnrm;
      i = 0;
      while ( i < curr_nev ) {	
	if (_ritzvaluesi[i] != zero) {
	  t_evecnrm = one/lapack.LAPY2(evecnrm[i],evecnrm[i+1]);
	  // Copy the real part of the eigenvector.  Scale by square-root of 2 to normalize the vector.
	  evecstempr = evecstemp->CloneView( index+i, 1 );
	  evecr1 = _evecs->CloneView( index+i, 1 );
	  evecr1->MvAddMv( t_evecnrm, *evecstempr, zero, *evecstempr );
	  delete evecr1; evecr1=0;
	  evecr1 = _evecs->CloneView( index+i+1, 1 );
	  evecr1->MvAddMv( t_evecnrm, *evecstempr, zero, *evecstempr );
	  // Note where imaginary part of eigenvector is.
	  indexi[conjprs] = i+1;
	  
	  // Increment counters.
	  conjprs++;
	  i = i+2; 				
	} else {
	  // Copy the real part of the eigenvector, scale to be norm one.
	  // We don't have to do anything for the imaginary
	  // part since we initialized the vectors to zero.
	  evecstempr = evecstemp->CloneView( index+i, 1 );
	  evecr1 = _evecs->CloneView( index+i, 1 );
	  evecr1->MvAddMv( one/evecnrm[i], *evecstempr, zero, *evecstempr );
	  // Increment counter.
	  i++;			
	}
	delete evecr1, evecr1=0;
	delete evecstempr; evecstempr=0;
      }
      // Set the imaginary part of the eigenvectors if conjugate pairs exist.
      // If the last eigenvector has a split conjugate pair, don't set negative imaginary
      // part.
      if (conjprs) {	
	MultiVec<TYPE>  *evecstempi=0, *eveci1=0;
	int* indexi_pnev = new int[ conjprs ];
	//
	// There is storage for an extra eigenvector.  
	// So, when the last eigenvalues is the first of a conjugate pair, that eigenvector will be computed.
	//
	for (i=0; i<conjprs; i++) {
	  indexi_pnev[i] = indexi[i] + _nev;
	  t_evecnrm = one/lapack.LAPY2(evecnrm[indexi[i]],evecnrm[indexi[i]-1]);
	  evecstempi = evecstemp->CloneView( indexi+i, 1 ); 
	  eveci1 = _evecs->CloneView( indexi_pnev+i, 1 );
	  eveci1->MvAddMv( t_evecnrm*Teuchos::ScalarTraits<TYPE>::magnitude(_ritzvaluesi[indexi[i]])/_ritzvaluesi[indexi[i]],
			   *evecstempi, zero, *evecstempi );
	  delete eveci1; eveci1=0;
	  // Change index and set non-conjugate part of imag eigenvector.
	  indexi[i]--;
	  indexi_pnev[i]--;
	  eveci1 = _evecs->CloneView( indexi_pnev+i, 1 );
	  eveci1->MvAddMv( t_evecnrm*Teuchos::ScalarTraits<TYPE>::magnitude(_ritzvaluesi[indexi[i]])/_ritzvaluesi[indexi[i]],
			   *evecstempi, zero, *evecstempi );
	  delete eveci1; eveci1=0;
	  delete evecstempi; evecstempi=0;
	}	      
	delete [] indexi_pnev;
      }
      
      // Clean up.
      delete evecstemp; 
      delete [] indexi; 
      delete [] evecnrm;
    }
    
    _isevecscurrent = true;
    delete [] index;
    delete basistemp;
  }
  
  template<class TYPE>
  void BlockKrylovSchur<TYPE>::ComputeSchurForm( const bool apply )
  {
    int m = _jstart*_block, n=_jstart*_block;
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    Teuchos::BLAS<int,TYPE> blas; 
    Teuchos::SerialDenseMatrix<int,TYPE> *Hj;
    Teuchos::SerialDenseMatrix<int,TYPE> Q( n, n );

    if (apply) {
      Hj = new Teuchos::SerialDenseMatrix<int,TYPE>( Teuchos::View, _hessmatrix, m, n );		
    } else {	
      Hj = new Teuchos::SerialDenseMatrix<int,TYPE>( Teuchos::Copy, _hessmatrix, m, n );
    }
    //
    SortSchurForm( *Hj, Q );
    //
    if (_nevblock <= _jstart || apply ) {
      //
      // Necessary variables.
      //
      int i=0,j=0       ;
      int mm1 = (_jstart-1)*_block;
      int _nevtemp, numimag;
      //
      // Determine new offset depending upon placement of conjugate pairs.	
      // ( if we are restarting, determine the new offset ) 
      if (apply) {
	_offset = _maxoffset;
	for (i=0; i<_maxoffset; i++) {
	  numimag = 0;
	  for (j=0; j<(_nevblock+i)*_block; j++) { 
	    if (_ritzvaluesi[j]!=zero) { numimag++; }; 
	  }
	  if (!(numimag % 2)) { _offset = i; break; }
	}
      }
      _nevtemp = n;
      if (_jstart > _nevblock+_offset)
	_nevtemp = (_nevblock+_offset)*_block;
      //
      Teuchos::SerialDenseMatrix<int,TYPE> sub_block_hess(Teuchos::View, _hessmatrix, _block, _block, m, mm1);
      Teuchos::SerialDenseMatrix<int,TYPE> sub_block_q(Teuchos::View, Q, _block, _nevtemp, mm1 );
      Teuchos::SerialDenseMatrix<int,TYPE> sub_block_b( _block, _nevtemp );
      blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, _block, _nevtemp, _block, one, 
		 sub_block_hess.values(), sub_block_hess.stride(), sub_block_q.values(), 
		 sub_block_q.stride(), zero, sub_block_b.values(), _block );
      //
      //---------------------------------------------------
      // Compute Schur decomposition error
      //
      // The residual for the Schur decomposition A(VQ) = (VQ)T + FB_m^TQ
      // where HQ = QT is || FB_m^TQ || = || H_{m+1,m}*B_m^TQ ||.
      //
      // We are only interested in the partial Krylov-Schur decomposition corresponding
      // to the _nev eigenvalues of interest or the _nevblock*_block number of
      // eigenvalues we're keeping.
      // NOTE:  The Schur error is not updated if the Schur decomposition is
      //        not large enough to compute _nev eigenvalues, else we could accidently
      //        satisfy a condition for convergence.
      //---------------------------------------------------
      //
      if (_nevblock <= _jstart ) {
	//
	Teuchos::SerialDenseMatrix<int,TYPE> sub_block_b2(Teuchos::View, sub_block_b, _block, _nev);
	_schurerror = sub_block_b2.normFrobenius();
	//
	// Determine whether we need to continue with the computations.
	//
	if (_schurerror < _residual_tolerance )
	  _exit_flg = true;
      }
      if (apply) {
	//
	//  We are going to restart, so update the Krylov-Schur decomposition.
	//
	// Update the Krylov-Schur basis.  
	//	
	int *index = new int[ n ]; assert(index!=NULL);
	for (i = 0; i < n; i++ ) {
	  index[i] = i;
	}
	Teuchos::SerialDenseMatrix<int,TYPE> Qnev(Teuchos::View, Q, n, _nevtemp);
	MultiVec<TYPE>* basistemp = _basisvecs->CloneView( index, _nevtemp );
	MultiVec<TYPE>* basistemp2 = _basisvecs->CloneCopy( index, n );
	basistemp->MvTimesMatAddMv ( one, *basistemp2, Qnev, zero );
	//
	// Update the Krylov-Schur quasi-triangular matrix.
	//
	Teuchos::SerialDenseMatrix<int,TYPE> Hjp1(Teuchos::View, _hessmatrix,_block,_nevtemp, _nevtemp );
	for (i=0; i<_block; i++) {
	  for (j=0; j<_nevtemp; j++) {
	    Hjp1(i, j) = sub_block_b(i, j);
	  }
	}      
	//
	// Clean up.
	//
	delete basistemp;
	delete basistemp2;
	delete [] index;
      }
    }
    //
    // Clean up.
    //
    delete Hj;
  }
  
  template<class TYPE>
  void BlockKrylovSchur<TYPE>::SortSchurForm( Teuchos::SerialDenseMatrix<int,TYPE>& H, Teuchos::SerialDenseMatrix<int,TYPE>& Q ) {
    const TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
    const TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    Teuchos::LAPACK<int,TYPE> lapack; 
    Teuchos::BLAS<int,TYPE> blas;
    int i, j, info=0;
    int n = H.numRows(), ldh = H.stride(), ldq = Q.stride(); 
    int m = H.numRows(), mm1 = H.numRows() - _block;
    TYPE* ptr_h = H.values();
    TYPE* ptr_q = Q.values();
    //
    //  If the operator is symmetric, analyze the block tridiagonal matrix
    //  and enforce symmetry.
    //
    if (_problem.IsSymmetric()) {
      if (_restartiter > 0 && _restarts!=0) {
	//
	// The method has been restarted, so more caution must be used in
	// imposing symmetry.
	//
	for(j=_nevblock*_block; j<n; j++) {
	  for(i=0; i<j; i++) {
	    H( i, j ) = H( j, i );
	  }
	}
      } else {
	//
	// We haven't restarted, so just enforce symmetry throughout Hj.
	//
	for( j=0; j<n; j++ ) {
	  for( i=0; i<j; i++ ) {
	    H( i, j ) = H( j, i );
	  }
	}
      }
    }
    //
    //---------------------------------------------------
    // Compute the current eigenvalue estimates
    // ---> Use driver GEES to first reduce to upper Hessenberg 
    // 	form and then compute Schur form, outputting eigenvalues
    //---------------------------------------------------
    //
    int lwork = 4*n;
    TYPE *work = new TYPE[lwork]; assert(work!=NULL);
    int *select = new int[ n ];
    int sdim = 0; 
    int *bwork = new int[ n ];
    char jobvs = 'V';
    char sort = 'N';
    lapack.GEES( jobvs, sort, select, n, ptr_h, ldh, &sdim,_ritzvalues,
		 _ritzvaluesi, ptr_q, ldq, work, lwork, bwork, &info );
    assert(info==0);
    //
    //---------------------------------------------------
    // Compute the current Ritz residuals for ALL the eigenvalues estimates (Ritz values)
    //           || Ax - x\theta || = || FB_m^Ts || 
    //                              = || V_m+1*H_{m+1,m}*B_m^T*s ||
    //                              = || H_{m+1,m}*B_m^T*s ||
    //
    // where V_m is the current Krylov-Schur basis and x = V_m*s
    // NOTE: This means that s = e_i if the problem is symmetric, else the eigenvectors
    //       of the Schur form need to be computed.
    //
    // First compute H_{m+1,m}*B_m^T, then determine what 's' is.
    //---------------------------------------------------
    //
    // H_{m+1,m}
    Teuchos::SerialDenseMatrix<int,TYPE> sub_block_hess(Teuchos::View, _hessmatrix, _block, _block, m, mm1);
    //
    // Last block rows of Q since the previous B_m is E_m (the last m-block of canonical basis vectors)
    Teuchos::SerialDenseMatrix<int,TYPE> sub_block_q(Teuchos::View, Q, _block, n, mm1 );
    //
    // Compute H_{m+1,m}*B_m^T
    Teuchos::SerialDenseMatrix<int,TYPE> sub_block_b( _block, n );
    blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, _block, n, _block, one, 
	       sub_block_hess.values(), sub_block_hess.stride(), sub_block_q.values(), 
	       sub_block_q.stride(), zero, sub_block_b.values(), _block );
    //
    // Determine what 's' is and compute Ritz residuals.
    //
    TYPE* b_ptr = sub_block_b.values();
    if (_problem.IsSymmetric()) {
      //
      // 's' is the i-th canonical basis vector.
      //
      for (i=0; i<n ; i++) {
	//_ritzresiduals[i] = blas.NRM2(_block, b_ptr + i*_block, 1);
	_ritzresiduals[i] = blas.NRM2(_block, b_ptr + i*_block, 1);
      }   
    } else {
      //
      //  's' is the eigenvector of the block upper triangular, Schur matrix.
      //
      char side = 'R';
      char howmny = 'A';
      int mm, ldvl = 1;
      TYPE *vl = new TYPE[ ldvl ];
      Teuchos::SerialDenseMatrix<int,TYPE> Q_temp( n, n );
      Teuchos::SerialDenseMatrix<int,TYPE> S( _block, n );
      lapack.TREVC( side, howmny, select, n, H.values(), H.stride(), vl, ldvl,
		  Q_temp.values(), Q_temp.stride(), n, &mm, work, &info );
      assert(info==0);
      delete [] vl;
      //
      // Scale the eigenvectors so that their euclidean norms are all one.
      // ( conjugate pairs get normalized by the sqrt(2) )
      //
      TYPE temp;
      TYPE* qt_ptr = Q_temp.values();
      i = 0;
      while( i < n ) {
	if ( _ritzvaluesi[i] != zero ) {
	  temp = lapack.LAPY2( blas.NRM2( n, qt_ptr+i*n, 1 ), blas.NRM2( n, qt_ptr+(i+1)*n, 1 ) );
	  blas.SCAL( n, one/temp, qt_ptr+i*n, 1 );
	  blas.SCAL( n, one/temp, qt_ptr+(i+1)*n, 1 );	      
	  i = i+2;
	} else {
	  temp = blas.NRM2( n, qt_ptr+i*n, 1 );
	  blas.SCAL( n, one/temp, qt_ptr+i*n, 1 );
	  i++;
	}
      }
      //
      // Compute H_{m+1,m}*B_m^T*S where the i-th column of S is 's' for the i-th Ritz-value
      //
      blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, _block, n, n, one, 
		 sub_block_b.values(), sub_block_b.stride(), Q_temp.values(), 
		 Q_temp.stride(), zero, S.values(), S.stride() );
      TYPE* s_ptr = S.values();
      i = 0;
      while( i < n ) {
	if ( _ritzvaluesi[i] != zero ) {
	  _ritzresiduals[i] = lapack.LAPY2( blas.NRM2(_block, s_ptr + i*_block, 1),
					    blas.NRM2(_block, s_ptr + (i+1)*_block, 1) );
	  _ritzresiduals[i+1] = _ritzresiduals[i];
	  i = i+2;
	} else {
	  _ritzresiduals[i] = blas.NRM2(_block, s_ptr + i*_block, 1);
	  i++;
	}
      }
    }
    //
    //---------------------------------------------------
    // Sort the eigenvalues
    //---------------------------------------------------
    //
    if (_problem.IsSymmetric())
      _sm.sort( this, n, _ritzvalues, _order );
    else
      _sm.sort( this, n, _ritzvalues, _ritzvaluesi, _order );
    //
    // Re-sort _ritzresiduals based on _order
    //
    TYPE* ritz2 = new TYPE[ n ];
    for (i=0; i<n; i++) { ritz2[i] = _ritzresiduals[ _order[i] ]; }
    blas.COPY( n, ritz2, 1, _ritzresiduals, 1 );
    delete [] ritz2;
    //
    // Copy the nev eigenvalues into the proper vectors
    // NOTE:  If we don't have nev Ritz values, then only n are copied
    //
    ( n > _nev ? blas.COPY( _nev, _ritzvalues, 1, _evals, 1 ) : blas.COPY( n, _ritzvalues, 1, _evals, 1 ) );
    if (!_problem.IsSymmetric() )
      ( n > _nev ? blas.COPY( _nev, _ritzvaluesi, 1, _evals+_nev, 1 ) : blas.COPY( n, _ritzvaluesi, 1, _evals+_nev, 1 ) );
    //
    //---------------------------------------------------
    // Reorder real Schur factorization, remember to add one to the indices for the
    // fortran call and determine offset.  The offset is necessary since the TREXC
    // method reorders in a nonsymmetric fashion, thus we use the reordering in
    // a stack-like fashion.  Also take into account conjugate pairs, which may mess
    // up the reordering, since the pair is moved if one of the pair is moved.
    //---------------------------------------------------
    //
    //cout<<"Before sorting the Schur form (H):"<<endl;
    //H.print(cout);	  
    //
    int _nevtemp = 0;
    char * compq = "V";
    int *offset2 = new int[ n ]; assert(offset2!=NULL);
    int *_order2 = new int[ n ]; assert(_order2!=NULL);
    i = 0; 
    while (i < n) {
      if (_ritzvaluesi[i] != zero) {
	offset2[_nevtemp] = 0;
	for (j=i; j<n; j++) {
	  if (_order[j] > _order[i]) { offset2[_nevtemp]++; }
	}
	_order2[_nevtemp] = _order[i];
	i = i+2;
      } else {
	offset2[_nevtemp] = 0;
	for (j=i; j<n; j++) {
	  if (_order[j] > _order[i]) { offset2[_nevtemp]++; }
	}
	_order2[_nevtemp] = _order[i];
	i++;
      }
      _nevtemp++;
    }
    for (i=_nevtemp-1; i>=0; i--) {
      lapack.TREXC( *compq, n, ptr_h, ldh, ptr_q, ldq, _order2[i]+1+offset2[i], 
		    1, work, &info );
      assert(info==0);
    }
    //cout<<"After sorting and reordering the Schur form (H):"<<endl;
    //H.print(cout); 
    //
    // Determine largest off diagonal element of Schur matrix for symmetric case.
    //
    TYPE _maxsymmelem = zero;
    if (_problem.IsSymmetric()) {
      for(j=0; j<n; j++){
	for(i=0; i<j; i++) {
	  if(Teuchos::ScalarTraits<TYPE>::magnitude(H(i, j))>_maxsymmelem) { _maxsymmelem = H(i, j); }
	}
      }
    }
    delete [] work; 
    delete [] bwork;
    delete [] select;
    delete [] offset2;
    delete [] _order2;
  }
  
  template<class TYPE>
  void BlockKrylovSchur<TYPE>::Restart() {
    //  This method assumes the ComputeResiduals has been called before it
    //  to compute the Schur vectors and residuals.  This information is used to 
    //  restart the factorization.
    //
    int i,j;
    int _nevtemp = (_nevblock+_offset)*_block;
    int *index = new int[ _nevtemp ];
    //
    //  Move the F_vec block to the _jstart+1 position.	
    //
    for (i=0; i<_block; i++) {
      index[i] = _jstart*_block + i;
    }
    MultiVec<TYPE>* F_vec = _basisvecs->CloneCopy(index, _block);
    for (i=0; i<_block; i++) {
      index[i] = _nevtemp + i;
    }
    _basisvecs->SetBlock( *F_vec, index, _block);
    //
    //  Reset the pointer.
    //
    _jstart = _nevblock+_offset; 
    //
    //  Clean up
    //
    delete F_vec;
    delete [] index;
  }
  
  template<class TYPE>
  void BlockKrylovSchur<TYPE>::CheckSchurVecs( const int j ) {
    //
    // Check the difference between the projection of A with the Schur vectors and the Schur matrix.
    // 
    int i, n = j*_block;
    int* index = new int[ n ];
    for( i=0; i<n; i++ ) { index[i] = i; } 
    TYPE one = Teuchos::ScalarTraits<TYPE>::one();
    Teuchos::SerialDenseMatrix<int,TYPE> Hj( Teuchos::View, _hessmatrix, n, n );
    Teuchos::SerialDenseMatrix<int,TYPE> SchurProj( n, n );
    MultiVec<TYPE>* Z = _basisvecs->CloneView( index, n );
    MultiVec<TYPE>* basistemp = _basisvecs->Clone( n );
    _Op->Apply( *Z, *basistemp );
    basistemp->MvTransMv( one, *Z, SchurProj );
    SchurProj.scale( -one );
    SchurProj += Hj;
    cout<< "Error in Schur Projection ( || (VQ)^T*A*(VQ) - S || ) at restart " << _restartiter << " is "<< SchurProj.normFrobenius()<<" (should be small)"<<endl;
  }
  
  
  template<class TYPE>
  void BlockKrylovSchur<TYPE>::CheckBlkArnRed( const int j ) {
    int i,k,m=(j+1)*_block;
    int *index = new int[m];
    ReturnType ret;       
    
    for ( i=0; i<_block; i++ ) {
      index[i] = m+i;
    }
    MultiVec<TYPE>* F_vec = _basisvecs->CloneView(index, _block);
    assert(F_vec!=NULL);
    
    TYPE *ptr_norms = new TYPE[m];
    TYPE sum=0.0;
    
    F_vec->MvNorm(ptr_norms);
    for ( i=0; i<_block; i++ ) {
      sum += ptr_norms[i];
    }
    
    for ( i=0; i<m; i++ ) {
      index[i] = i;  
    }
    MultiVec<TYPE>* Vj = _basisvecs->CloneView(index, m);
    assert(Vj!=NULL);   
    cout << " " << endl;
    cout << "********Block Arnoldi iteration******** " << j+1 << endl;
    cout << " " << endl;
    
    const TYPE one=1.0;
    const TYPE zero=0.0;
    Teuchos::SerialDenseMatrix<int,TYPE> VTV(m,m);
    ret = _problem.InnerProd( *Vj, *Vj, VTV );
    if (ret != Ok) { }
    TYPE* ptr=VTV.values();
    TYPE column_sum;
    
    for (k=0; k<m; k++) {
      column_sum=zero;
      for (i=0; i<m; i++) {
	if (i==k) {
	  ptr[i] -= one;
	}
	column_sum += ptr[i];
      }
      cout <<  " V^T*B*V-I " << "for column " << k << " is " << Teuchos::ScalarTraits<TYPE>::magnitude(column_sum) << endl;
      ptr += m;
    }
    cout << " " << endl;
    
    Teuchos::SerialDenseMatrix<int,TYPE> E(m,_block);
    
    ret = _problem.InnerProd( *Vj, *F_vec, E );
    if (ret != Ok) { }
    TYPE* ptr_Ej=E.values();
    
    for (k=0;k<_block;k++) {
      column_sum=zero;
      for (i=0; i<m; i++) {
	column_sum += ptr_Ej[i];
      }
      ptr_Ej += m;
      if (ptr_norms[k]) column_sum = column_sum/ptr_norms[k];
      cout << " B-Orthogonality with F " << "for column " << k << " is " << Teuchos::ScalarTraits<TYPE>::magnitude(column_sum) << endl;
    }
    cout << " " << endl;
                 
    MultiVec<TYPE>* AVj = _basisvecs->Clone(m); assert(AVj!=NULL);
    ret = _Op->Apply(*Vj,*AVj);
    Teuchos::SerialDenseMatrix<int,TYPE> Hj(Teuchos::View, _hessmatrix, m, m);
    AVj->MvTimesMatAddMv(-one, *Vj, Hj, one);
    for ( i=0; i<_block; i++ ) {  
      index[i] = j*_block+i;
    }
    
    MultiVec<TYPE>* Fj = AVj->CloneView(index, _block);
    Fj->MvAddMv(-one, *F_vec, one, *Fj);
    
    AVj->MvNorm(ptr_norms);
    
    for ( i=0; i<m; i++ ) { 
      cout << " Arnoldi relation " << "for column " << i << " is " << ptr_norms[i] << endl;  
    }
    cout << " " << endl;
    
    delete F_vec;
    delete Fj;
    delete AVj;
    delete Vj;
    delete [] index;
    delete [] ptr_norms;
  }
  
} // End of namespace Anasazi
#endif
// End of file AnasaziBlockKrylovSchur.hpp



