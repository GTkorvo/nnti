// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef QP_SOLVER_RELAXED_QPOPTSOL_H
#define QP_SOLVER_RELAXED_QPOPTSOL_H

#include <vector>

#include "ConstrainedOptPack_QPSolverRelaxed.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Node base clase for the primal QP solvers QPOPT and QPSOL.
 *
 * In this implementation it is required that <tt>G</tt> only support
 * the <tt>MatrixOp</tt> interface and is therefore quite flexible in
 * the QPs it can solve.
 */
class QPSolverRelaxedQPOPTSOL : public QPSolverRelaxed
{
public:

  // /////////////////////////////////////
  /** @name Public Types */
  //@{

  /** \brief . */
  typedef FortranTypes::f_int       f_int;
  /** \brief . */
  typedef FortranTypes::f_dbl_prec  f_dbl_prec;
  /** \brief . */
  typedef FortranTypes::f_logical   f_logical;

  //@}

  /** \brief . */
  QPSolverRelaxedQPOPTSOL();

  /** \brief . */
  ~QPSolverRelaxedQPOPTSOL();

  /// Return a pointer to the matrix G to be used in the calculation of H*x by QPOPT and QPSOL.
  virtual const MatrixOp* G() const;

  /// Return the value of the "big M" used in the relaxation (called by QPHESS functions).
  virtual value_type use_as_bigM() const;

  // /////////////////////////////////
  // Overridden from QPSolverRelaxed

  /** \brief . */
  QPSolverStats get_qp_stats() const;

  /** \brief . */
  void release_memory();

protected:

  // /////////////////////////////////
  // Overridden from QPSolverRelaxed

  /** \brief . */
  QPSolverStats::ESolutionType imp_solve_qp(
    std::ostream* out, EOutputLevel olevel, ERunTests test_what
    ,const Vector& g, const MatrixSymOp& G
    ,value_type etaL
    ,const Vector* dL, const Vector* dU
    ,const MatrixOp* E, BLAS_Cpp::Transp trans_E, const Vector* b
    ,const Vector* eL, const Vector* eU
    ,const MatrixOp* F, BLAS_Cpp::Transp trans_F, const Vector* f
    ,value_type* obj_d
    ,value_type* eta, VectorMutable* d
    ,VectorMutable* nu
    ,VectorMutable* mu, VectorMutable* Ed
    ,VectorMutable* lambda, VectorMutable* Fd
    );

  // //////////////////////////////////////////////////////////////
  // Protected types

  /** \brief . */
  typedef std::vector<f_int>			ISTATE_t;
  /** \brief . */
  typedef std::vector<f_int>			IWORK_t;
  /** \brief . */
  typedef std::vector<f_dbl_prec>		WORK_t;
public: // RAB: 2001/05/03: MS VC++ 6.0 must have this public ???
  /** \brief . */
  enum EInform {
    STRONG_LOCAL_MIN,
    WEAK_LOCAL_MIN,
    MAX_ITER_EXCEEDED,
    OTHER_ERROR
  };

protected:

  // //////////////////////////////////////////////////////////////
  // Protected Data Members.

  QPSolverStats			qp_stats_;

  /** @name Input/Output parameters common to both QPOPT and QPSOL
    *
    * These are access and updated by subclasses that call QPOPT
    * and QPSOL.
    */
  //@{

  /** \brief . */
  f_int		N_;
  /** \brief . */
  f_int		NCLIN_;
  /** \brief . */
  DMatrix	A_;
  /** \brief . */
  DVector		BL_;
  /** \brief . */
  DVector		BU_;
  /** \brief . */
  DVector		CVEC_;
  /** \brief . */
  ISTATE_t	ISTATE_;
  /** \brief . */
  DVector		X_;
  /** \brief . */
  DVector		AX_;
  /** \brief . */
  DVector		CLAMDA_;
  /** \brief . */
  f_int		ITER_;
  /** \brief . */
  f_dbl_prec	OBJ_;
  /** \brief . */
  f_int		LIWORK_;
  /** \brief . */
  IWORK_t		IWORK_;
  /** \brief . */
  f_int		LWORK_;
  /** \brief . */
  WORK_t		WORK_;

  //@}

  // /////////////////////////////////////////////////////////////
  // Template method primatives to be overridden.

  /// Length of integer workspace
  virtual f_int liwork(f_int N, f_int NCLIN) const = 0;

  /// Length of real workspace
  virtual f_int lrwork(f_int N, f_int NCLIN) const = 0;

  /** \brief Solve the QP defined in the protected input data members
    * and set the solution in the protected output data members.
    */
  virtual EInform call_qp_solver(bool warm_start) = 0;

private:

  // /////////////////////////
  // Private types

  typedef std::vector<f_int>	ibnds_t;

  // ///////////////////////////
  // Private data members

  size_type			n_inequ_bnds_;		// Used to record the number of bounds with at least
                      // one bound existing in eL and eU.
  ibnds_t				i_inequ_bnds_;		// size(nbounds_). Remembers which bounds in
                      // eL, eU had at least one bound present.  This is
                      // needed to map from CLAMDA_ to mu.
  value_type			bigM_;				// Big M value used to construct relaxation.
  value_type			use_as_bigM_;		// Big M value used in QPHESS.
  const MatrixOp*	G_;					// used to compute HESS * x = [ G, 0; 0, bigM ] * x products.

  // ///////////////////////////
  // Private member functions

  // not defined and not to be called.
  QPSolverRelaxedQPOPTSOL(const QPSolverRelaxedQPOPTSOL&);
  QPSolverRelaxedQPOPTSOL& operator=(const QPSolverRelaxedQPOPTSOL&);

};	// end class QPSolverRelaxedQPOPTSOL

}	// end namespace ConstrainedOptimizationPackTypes

#endif	// QP_SOLVER_RELAXED_QPOPTSOL_H
