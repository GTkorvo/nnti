// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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

#ifndef EPETRA_EXT_MODEL_EVALUATOR_HPP
#define EPETRA_EXT_MODEL_EVALUATOR_HPP

#include "EpetraExt_ConfigDefs.h"
#include "EpetraExt_PolynomialVectorTraits.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_Polynomial.hpp"
#include "Teuchos_Array.hpp"

class Epetra_Map;
class Epetra_Vector;
class Epetra_Operator;

namespace EpetraExt {

/** \brief Base interface for evaluating a stateless "model".
 *
 * ToDo: Finish Documentation!
 */
class ModelEvaluator : virtual public Teuchos::Describable {
public:

  /** \name Public types */
  //@{

  /** \brief.  */
  enum EInArgsMembers {
    IN_ARG_x_dot
    ,IN_ARG_x
    ,IN_ARG_x_dot_poly 
    ,IN_ARG_x_poly 
    ,IN_ARG_t
    ,IN_ARG_alpha
    ,IN_ARG_beta
  };
  static const int NUM_E_IN_ARGS_MEMBERS=7;

  /** \brief . */
  class InArgs {
  public:
    /** \brief. */
    InArgs();
    /** \brief . */
    std::string modelEvalDescription() const;
    /** \brief .  */
    int Np() const;
    /** \brief. */
    void set_x_dot( const Teuchos::RefCountPtr<const Epetra_Vector> &x_dot );
    /** \brief. */
    Teuchos::RefCountPtr<const Epetra_Vector> get_x_dot() const;
    /** \brief. */
    void set_x( const Teuchos::RefCountPtr<const Epetra_Vector> &x );
    /** \brief. */
    Teuchos::RefCountPtr<const Epetra_Vector> get_x() const;
    void set_x_poly(
      const Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > &x_poly
      );
    /** \brief .  */
    Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > get_x_poly() const;
    /** \brief .  */
    void set_x_dot_poly(
      const Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > &x_dot_poly
      );
    /** \brief .  */
    Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > get_x_dot_poly() const;
    /** \brief. */
    void set_p( int l, const Teuchos::RefCountPtr<const Epetra_Vector> &p_l );
    /** \brief. */
    Teuchos::RefCountPtr<const Epetra_Vector> get_p(int l) const;
    /** \brief. */
    void set_t( double t );
    /** \brief. */
    double get_alpha() const;
    /** \brief. */
    void set_alpha( double alpha );
    /** \brief. */
    double get_beta() const;
    /** \brief. */
    void set_beta( double beta );
    /** \brief. */
    double get_t() const;
    /** \brief. */
    bool supports(EInArgsMembers arg) const;
  protected:
    /** \brief . */
    void _setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void _set_Np(int Np);
    /** \brief . */
    void _setSupports( EInArgsMembers arg, bool supports );
  private:
    // types
    typedef Teuchos::Array<Teuchos::RefCountPtr<const Epetra_Vector> > p_t;
    // data
    std::string modelEvalDescription_;
    Teuchos::RefCountPtr<const Epetra_Vector>  x_dot_;
    Teuchos::RefCountPtr<const Epetra_Vector>  x_;
    Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > x_dot_poly_;
    Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > x_poly_;
    p_t                                        p_;
    double                                     t_;
    double                                     alpha_;
    double                                     beta_;
    bool supports_[NUM_E_IN_ARGS_MEMBERS];
    // functions
    void assert_supports(EInArgsMembers arg) const;
    void assert_l(int l) const;
  };

  /** \brief. */
  enum EEvalType {
    EVAL_TYPE_EXACT                ///< Exact function evaluation
    ,EVAL_TYPE_APPROX_DERIV        ///< An approximate derivative (i.e. for a Jacobian)
    ,EVAL_TYPE_VERY_APPROX_DERIV   ///< A very approximate derivative (i.e. for a preconditioner)
  };

  /** \brief . */
  template<class ObjType>
  class Evaluation : public Teuchos::RefCountPtr<ObjType> {
  public:
    /** \brief . */
    Evaluation() : evalType_(EVAL_TYPE_EXACT) {}
    /** \brief . */
    Evaluation( const Teuchos::RefCountPtr<ObjType> &obj )
      : Teuchos::RefCountPtr<ObjType>(obj), evalType_(EVAL_TYPE_EXACT) {}
    /** \brief . */
    Evaluation( const Teuchos::RefCountPtr<ObjType> &obj, EEvalType evalType )
      : Teuchos::RefCountPtr<ObjType>(obj), evalType_(evalType) {}
    /** \brief . */
    EEvalType getType() const { return evalType_; }
    /** \brief . */
    void reset( const Teuchos::RefCountPtr<ObjType> &obj, EEvalType evalType ) 
    { this->operator=(obj); evalType_ = evalType; }
  private:
    EEvalType                      evalType_;
  };
  
  /** \brief . */
  enum EDerivativeMultiVectorOrientation {
    DERIV_MV_BY_COL           ///< .
    ,DERIV_TRANS_MV_BY_ROW    ///< .
  };

  /** \brief . */
  enum EDerivativeLinearOp { DERIV_LINEAR_OP };

  /** \brief . */
  class DerivativeSupport {
  public:
    /** \brief . */
    DerivativeSupport()
      :supportsLinearOp_(false), supportsMVByCol_(false), supportsTransMVByRow_(false)
      {}
    /** \brief . */
    DerivativeSupport( EDerivativeLinearOp )
      :supportsLinearOp_(true), supportsMVByCol_(false), supportsTransMVByRow_(false)
      {}
    /** \brief . */
    DerivativeSupport( EDerivativeMultiVectorOrientation mvOrientation )
      :supportsLinearOp_(false), supportsMVByCol_(mvOrientation==DERIV_MV_BY_COL)
      ,supportsTransMVByRow_(mvOrientation==DERIV_TRANS_MV_BY_ROW)
      {}
    /** \brief . */
    DerivativeSupport(
      EDerivativeLinearOp, EDerivativeMultiVectorOrientation mvOrientation )
      :supportsLinearOp_(true), supportsMVByCol_(mvOrientation==DERIV_MV_BY_COL)
      ,supportsTransMVByRow_(mvOrientation==DERIV_TRANS_MV_BY_ROW)
      {}
    /** \brief . */
    DerivativeSupport(
      EDerivativeMultiVectorOrientation mvOrientation1,
      EDerivativeMultiVectorOrientation mvOrientation2
      )
      :supportsLinearOp_(false)
      ,supportsMVByCol_(
        mvOrientation1==DERIV_MV_BY_COL||mvOrientation2==DERIV_MV_BY_COL )
      ,supportsTransMVByRow_(
        mvOrientation1==DERIV_TRANS_MV_BY_ROW||mvOrientation2==DERIV_TRANS_MV_BY_ROW )
      {}
    /** \brief . */
    DerivativeSupport& plus(EDerivativeLinearOp)
      { supportsLinearOp_ = true; return *this; }
    /** \brief . */
    DerivativeSupport& plus(EDerivativeMultiVectorOrientation mvOrientation)
      {
        switch(mvOrientation) {
          case DERIV_MV_BY_COL: supportsMVByCol_ = true; break;
          case DERIV_TRANS_MV_BY_ROW: supportsTransMVByRow_ = true; break;
          default: TEST_FOR_EXCEPT(true);
        }
        return *this;
      }
    /** \brief . */
    bool none() const
      { return ( !supportsLinearOp_ && !supportsMVByCol_ && !supportsTransMVByRow_ ); }
    /** \brief . */
    bool supports(EDerivativeLinearOp) const
      { return supportsLinearOp_; }
    /** \brief . */
    bool supports(EDerivativeMultiVectorOrientation mvOrientation) const
      {
        switch(mvOrientation) {
          case DERIV_MV_BY_COL: return supportsMVByCol_;
          case DERIV_TRANS_MV_BY_ROW: return supportsTransMVByRow_;
          default: TEST_FOR_EXCEPT(true);
        }
        return false; // Will never be called!
      }
  private:
    bool supportsLinearOp_;
    bool supportsMVByCol_;
    bool supportsTransMVByRow_;
  public:
  };

  /** \brief . */
  enum EDerivativeLinearity {
    DERIV_LINEARITY_UNKNOWN      ///< .
    ,DERIV_LINEARITY_CONST       ///< .
    ,DERIV_LINEARITY_NONCONST    ///< .
  };
  /** \brief . */
  enum ERankStatus {
    DERIV_RANK_UNKNOWN       ///< .
    ,DERIV_RANK_FULL         ///< .
    ,DERIV_RANK_DEFICIENT    ///< .
  };

  /** \brief . */
  struct DerivativeProperties {
    /** \brief . */
    EDerivativeLinearity linearity;
    /** \brief . */
    ERankStatus rank;
    /** \brief . */
    bool supportsAdjoint;
    /** \brief . */
    DerivativeProperties()
      :linearity(DERIV_LINEARITY_UNKNOWN),rank(DERIV_RANK_UNKNOWN),supportsAdjoint(false) {}
    /** \brief . */
    DerivativeProperties(
      EDerivativeLinearity in_linearity, ERankStatus in_rank, bool in_supportsAdjoint
      ):linearity(in_linearity),rank(in_rank),supportsAdjoint(in_supportsAdjoint) {}
  };

  /** \brief Simple aggregate class for a derivative object represented as a
   * column-wise multi-vector or its transpose as a row-wise multi-vector.
   */
  class DerivativeMultiVector {
  public:
    /** \brief . */
    DerivativeMultiVector() {}
    /** \brief . */
    DerivativeMultiVector(
      const Teuchos::RefCountPtr<Epetra_MultiVector> &mv
      ,const EDerivativeMultiVectorOrientation orientation = DERIV_MV_BY_COL
      ,const Teuchos::Array<int> &paramIndexes = Teuchos::Array<int>()
      ) : mv_(mv), orientation_(orientation), paramIndexes_(paramIndexes) {}
    /** \brief . */
    void changeOrientation( const EDerivativeMultiVectorOrientation orientation )
      { orientation_ = orientation; };
    /** \brief . */
    Teuchos::RefCountPtr<Epetra_MultiVector> getMultiVector() const
      { return mv_; }
    /** \brief . */
    EDerivativeMultiVectorOrientation getOrientation() const
      { return orientation_; }
    /** \brief . */
    const Teuchos::Array<int>& getParamIndexes() const
      { return paramIndexes_; }
  private:
    Teuchos::RefCountPtr<Epetra_MultiVector> mv_;
    EDerivativeMultiVectorOrientation orientation_;
    Teuchos::Array<int> paramIndexes_;
  };

  /** \brief Simple aggregate class that stores a derivative object
   * as a general linear operator or as a multi-vector.
   */
  class Derivative {
  public:
    /** \brief . */
    Derivative() {}
    /** \brief . */
    Derivative( const Teuchos::RefCountPtr<Epetra_Operator> &lo )
      : lo_(lo) {}
     /** \brief . */
    Derivative(
      const Teuchos::RefCountPtr<Epetra_MultiVector> &mv
      ,const EDerivativeMultiVectorOrientation orientation = DERIV_MV_BY_COL
      ) : dmv_(mv,orientation) {}
   /** \brief . */
    Derivative( const DerivativeMultiVector &dmv )
      : dmv_(dmv) {}
    /** \brief . */
    Teuchos::RefCountPtr<Epetra_Operator> getLinearOp() const
      { return lo_; }
    /** \brief . */
    Teuchos::RefCountPtr<Epetra_MultiVector> getMultiVector() const
      { return dmv_.getMultiVector(); }
    /** \brief . */
    EDerivativeMultiVectorOrientation getMultiVectorOrientation() const
      { return dmv_.getOrientation(); }
    /** \brief . */
    DerivativeMultiVector getDerivativeMultiVector() const
      { return dmv_; }
    /** \brief . */
    bool isEmpty() const
        { return !lo_.get() && !dmv_.getMultiVector().get(); }
  private:
    Teuchos::RefCountPtr<Epetra_Operator> lo_;
    DerivativeMultiVector dmv_;
  };

  /** \brief.  */
  enum EOutArgsMembers {
    OUT_ARG_f
    ,OUT_ARG_W
    ,OUT_ARG_f_poly
  };
  static const int NUM_E_OUT_ARGS_MEMBERS=3;

  /** \brief . */
  enum EOutArgsDfDp {
    OUT_ARG_DfDp   ///< .
  };

  /** \brief . */
  enum EOutArgsDgDx {
    OUT_ARG_DgDx   ///< .
  };

  /** \brief . */
  enum EOutArgsDgDp {
    OUT_ARG_DgDp   ///< .
  };

  /** \brief . */
  class OutArgs {
  public:
    /** \brief. */
    OutArgs();
    /** \brief . */
    std::string modelEvalDescription() const;
    /** \brief .  */
    int Np() const;
    /** \brief .  */
    int Ng() const;
    /** \brief. */
    bool supports(EOutArgsMembers arg) const;
    /** \brief <tt>0 <= l && l < Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDfDp arg, int l) const;
    /** \brief <tt>0 <= j && j < Ng()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDx arg, int j) const;
    /** \brief <tt>0 <= j && j < Ng()</tt> and <tt>0 <= l && l < Np()</tt>.  */
    const DerivativeSupport& supports(EOutArgsDgDp arg, int j, int l) const;
    /** \brief. */
    void set_f( const Evaluation<Epetra_Vector> &f );
    /** \brief. */
    Evaluation<Epetra_Vector> get_f() const;
    /** \brief Set <tt>g(j)</tt> where <tt>0 <= j && j < this->Ng()</tt>.  */
    void set_g( int j, const Evaluation<Epetra_Vector> &g_j );
    /** \brief Get <tt>g(j)</tt> where <tt>0 <= j && j < this->Ng()</tt>.  */
    Evaluation<Epetra_Vector> get_g(int j) const;
    /** \brief. */
    void set_W( const Teuchos::RefCountPtr<Epetra_Operator> &W );
    /** \brief. */
    Teuchos::RefCountPtr<Epetra_Operator> get_W() const;
    /** \brief . */
    DerivativeProperties get_W_properties() const;
    /** \brief .  */
    void set_DfDp(int l,  const Derivative &DfDp_l);
    /** \brief .  */
    Derivative get_DfDp(int l) const;
    /** \brief . */
    DerivativeProperties get_DfDp_properties(int l) const;
    /** \brief .  */
    void set_DgDx(int j, const Derivative &DgDx_j);
    /** \brief .  */
    Derivative get_DgDx(int j) const;
    /** \brief . */
    DerivativeProperties get_DgDx_properties(int j) const;
    /** \brief .  */
    void set_DgDp( int j, int l, const Derivative &DgDp_j_l );
    /** \brief .  */
    Derivative get_DgDp(int j, int l) const;
    /** \brief . */
    DerivativeProperties get_DgDp_properties(int j, int l) const;
    /** \brief .  */
    void set_f_poly( const Teuchos::RefCountPtr<Teuchos::Polynomial<Epetra_Vector> > &f_poly );
    /** \brief .  */
    Teuchos::RefCountPtr<Teuchos::Polynomial<Epetra_Vector> > get_f_poly() const;
    /** \brief Return true if the function or its derivatives are set. */
    bool funcOrDerivesAreSet(EOutArgsMembers arg) const;
  protected:
    /** \brief . */
    void _setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void _set_Np_Ng(int Np, int Ng);
    /** \brief . */
    void _setSupports( EOutArgsMembers arg, bool supports );
    /** \brief . */
    void _setSupports( EOutArgsDfDp arg, int l, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDx arg, int j, const DerivativeSupport& );
    /** \brief . */
    void _setSupports( EOutArgsDgDp arg, int j, int l, const DerivativeSupport& );
    /** \brief . */
    void _set_W_properties( const DerivativeProperties &W_properties );
    /** \brief . */
    void _set_DfDp_properties( int l, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDx_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void _set_DgDp_properties( int j, int l, const DerivativeProperties &properties );
  private:
    // types
    typedef Teuchos::Array<Evaluation<Epetra_Vector> > g_t;
    typedef Teuchos::Array<Derivative> deriv_t;
    typedef Teuchos::Array<DerivativeProperties> deriv_properties_t;
    typedef Teuchos::Array<DerivativeSupport> supports_t;
    // data
    std::string modelEvalDescription_;
    bool supports_[NUM_E_OUT_ARGS_MEMBERS];
    supports_t supports_DfDp_; // Np
    supports_t supports_DgDx_; // Ng
    supports_t supports_DgDp_; // Ng x Np
    Evaluation<Epetra_Vector> f_;
    g_t g_;
    Teuchos::RefCountPtr<Epetra_Operator> W_;
    DerivativeProperties W_properties_;
    deriv_t DfDp_; // Np
    deriv_properties_t DfDp_properties_; // Np
    deriv_t DgDx_; // Ng
    deriv_properties_t DgDx_properties_; // Ng
    deriv_t DgDp_; // Ng x Np
    deriv_properties_t DgDp_properties_; // Ng x Np
    Teuchos::RefCountPtr<Teuchos::Polynomial<Epetra_Vector> > f_poly_;
    // functions
    void assert_supports(EOutArgsMembers arg) const;
    void assert_supports(EOutArgsDfDp arg, int l) const;
    void assert_supports(EOutArgsDgDx arg, int j) const;
    void assert_supports(EOutArgsDgDp arg, int j, int l) const;
    void assert_l(int l) const;
    void assert_j(int j) const;
  };

  //@}

  /** \name Destructor */
  //@{

  /** \brief . */
  virtual ~ModelEvaluator();

  //@}

  /** \name Vector maps */
  //@{

  /** \breif . */
  virtual Teuchos::RefCountPtr<const Epetra_Map> get_x_map() const = 0;

  /** \breif . */
  virtual Teuchos::RefCountPtr<const Epetra_Map> get_f_map() const = 0;

  /** \breif . */
  virtual Teuchos::RefCountPtr<const Epetra_Map> get_p_map(int l) const;

  /** \brief Get the names of the parameters associated with parameter
   * subvector l if available.
   *
   * \return Returns an RCP to a Teuchos::Array<std::string> object that
   * contains the names of the parameters.  If returnVal == Teuchos::null,
   * then there are no names available for the parameter subvector p(l).  If
   * returnVal->size() == 1, then a single name is given to the entire
   * parameter subvector.  If returnVal->size() ==
   * get_p_map(l)->GlobalNumElements(), then a name is given to every
   * parameter scalar entry.
   *
   * The default implementation return returnVal==Teuchos::null which means
   * by default, parameters have no names associated with them.
   */
  virtual Teuchos::RefCountPtr<const Teuchos::Array<std::string> > get_p_names(int l) const;

  /** \breif . */
  virtual Teuchos::RefCountPtr<const Epetra_Map> get_g_map(int j) const;

  //@}

  /** \name Initial guesses for variables/parameters */
  //@{

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_x_init() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_x_dot_init() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_p_init(int l) const;

  /** \brief . */
  virtual double get_t_init() const;

  //@}

  /** \name Bounds for variables/parameters */
  //@{

  /** \brief Return the value of an infinite bound.
   *
   * The default implementation returns 1e+50.
   */
  virtual double getInfBound() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_x_lower_bounds() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_x_upper_bounds() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_p_lower_bounds(int l) const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<const Epetra_Vector> get_p_upper_bounds(int l) const;

  /** \brief . */
  virtual double get_t_lower_bound() const;

  /** \brief . */
  virtual double get_t_upper_bound() const;

  //@}

  /** \name Factory functions for creating derivative objects */
  //@{

  /** \brief If supported, create a <tt>Epetra_Operator</tt> object for
   * <tt>W</tt> to be evaluated.
   *
   * The default implementation returns <tt>return.get()==NULL</tt>
   * (i.e. implicit solvers are not supported by default).
   */
  virtual Teuchos::RefCountPtr<Epetra_Operator> create_W() const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<Epetra_Operator> create_DfDp_op(int l) const;

  // ToDo: Add functions for creating D(g(j))/D(x_dot) if needed!

  /** \brief . */
  virtual Teuchos::RefCountPtr<Epetra_Operator> create_DgDx_op(int j) const;

  /** \brief . */
  virtual Teuchos::RefCountPtr<Epetra_Operator> create_DgDp_op( int j, int l ) const;

  //@}

  /** \name Computational functions */
  //@{

  /** \brief . */
  virtual InArgs createInArgs() const = 0;

  /** \brief . */
  virtual OutArgs createOutArgs() const = 0;

  /** \brief . */
  virtual void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const = 0;

  //@}

protected:

  /** \name Protected types */
  //@{

  /** \brief . */
  class InArgsSetup : public InArgs {
  public:
    /** \brief . */
    void setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void set_Np(int Np);
    /** \brief . */
    void setSupports( EInArgsMembers arg, bool supports = true );
  };

  /** \brief . */
  class OutArgsSetup : public OutArgs {
  public:
    /** \brief . */
    void setModelEvalDescription( const std::string &modelEvalDescription );
    /** \brief . */
    void set_Np_Ng(int Np, int Ng);
    /** \brief . */
    void setSupports( EOutArgsMembers arg, bool supports = true );
    /** \brief . */
    void setSupports(EOutArgsDfDp arg, int l, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDx arg, int j, const DerivativeSupport& );
    /** \brief . */
    void setSupports(EOutArgsDgDp arg, int j, int l, const DerivativeSupport& );
    /** \brief . */
    void set_W_properties( const DerivativeProperties &properties );
    /** \brief . */
    void set_DfDp_properties( int l, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDx_properties( int j, const DerivativeProperties &properties );
    /** \brief . */
    void set_DgDp_properties( int j, int l, const DerivativeProperties &properties );
  };

  //@}

};

// ////////////////////////////
// Helper functions

/** \brief . */
std::string toString( ModelEvaluator::EDerivativeMultiVectorOrientation orientation );

/** \brief . */
std::string toString( ModelEvaluator::EInArgsMembers inArg );

/** \brief . */
std::string toString( ModelEvaluator::EOutArgsMembers outArg );

/** \brief . */
Teuchos::RefCountPtr<Epetra_Operator>
getLinearOp(
  const std::string &modelEvalDescription,
  const ModelEvaluator::Derivative &deriv,
  const std::string &derivName
  );

/** \brief . */
Teuchos::RefCountPtr<Epetra_MultiVector>
getMultiVector(
  const std::string &modelEvalDescription,
  const ModelEvaluator::Derivative &deriv,
  const std::string &derivName,
  const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  );

/** \brief . */
Teuchos::RefCountPtr<Epetra_Operator>
get_DfDp_op(
  const int l
  ,const ModelEvaluator::OutArgs &outArgs
  );

/** \brief . */
Teuchos::RefCountPtr<Epetra_MultiVector>
get_DfDp_mv(
  const int l
  ,const ModelEvaluator::OutArgs &outArgs
  );

/** \brief . */
Teuchos::RefCountPtr<Epetra_MultiVector>
get_DgDx_mv(
  const int j
  ,const ModelEvaluator::OutArgs &outArgs
  ,const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  );

/** \brief . */
Teuchos::RefCountPtr<Epetra_MultiVector>
get_DgDp_mv(
  const int j
  ,const int l
  ,const ModelEvaluator::OutArgs &outArgs
  ,const ModelEvaluator::EDerivativeMultiVectorOrientation mvOrientation
  );

// ///////////////////////////
// Inline Functions

//
// ModelEvaluator::InArgs
//

inline
std::string ModelEvaluator::InArgs::modelEvalDescription() const
{ return modelEvalDescription_; }

inline
int ModelEvaluator::InArgs::Np() const
{ return p_.size(); }

inline
void ModelEvaluator::InArgs::set_x_dot( const Teuchos::RefCountPtr<const Epetra_Vector> &x_dot )
{ assert_supports(IN_ARG_x_dot); x_dot_ = x_dot; }

inline
Teuchos::RefCountPtr<const Epetra_Vector> ModelEvaluator::InArgs::get_x_dot() const
{ assert_supports(IN_ARG_x_dot); return x_dot_; }

inline
void ModelEvaluator::InArgs::set_x( const Teuchos::RefCountPtr<const Epetra_Vector> &x )
{ assert_supports(IN_ARG_x); x_ = x; }

inline
Teuchos::RefCountPtr<const Epetra_Vector> ModelEvaluator::InArgs::get_x() const
{ assert_supports(IN_ARG_x); return x_; }

inline 
void ModelEvaluator::InArgs::set_x_dot_poly( const Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > &x_dot_poly )
{ assert_supports(IN_ARG_x_dot_poly); x_dot_poly_ = x_dot_poly; }

inline 
Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> >
ModelEvaluator::InArgs::get_x_dot_poly() const
{ assert_supports(IN_ARG_x_dot_poly); return x_dot_poly_; }

inline 
void ModelEvaluator::InArgs::set_x_poly( const Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> > &x_poly )
{ assert_supports(IN_ARG_x_poly); x_poly_ = x_poly; }

inline 
Teuchos::RefCountPtr<const Teuchos::Polynomial<Epetra_Vector> >
ModelEvaluator::InArgs::get_x_poly() const
{ assert_supports(IN_ARG_x_poly); return x_poly_; }

inline
void ModelEvaluator::InArgs::set_p( int l, const Teuchos::RefCountPtr<const Epetra_Vector> &p_l )
{ assert_l(l); p_[l] = p_l; }

inline
Teuchos::RefCountPtr<const Epetra_Vector> ModelEvaluator::InArgs::get_p(int l) const
{ assert_l(l); return p_[l]; }

inline
void ModelEvaluator::InArgs::set_t( double t )
{ assert_supports(IN_ARG_t); t_ = t; }

inline
double ModelEvaluator::InArgs::get_t() const
{ assert_supports(IN_ARG_t); return t_; }

inline
void ModelEvaluator::InArgs::set_alpha( double alpha )
{ assert_supports(IN_ARG_alpha); alpha_ = alpha; }

inline
double ModelEvaluator::InArgs::get_alpha() const
{ assert_supports(IN_ARG_alpha); return alpha_; }

inline
void ModelEvaluator::InArgs::set_beta( double beta )
{ assert_supports(IN_ARG_beta); beta_ = beta; }

inline
double ModelEvaluator::InArgs::get_beta() const
{ assert_supports(IN_ARG_beta); return beta_; }

inline
void ModelEvaluator::InArgs::_setModelEvalDescription( const std::string &modelEvalDescription )
{
  modelEvalDescription_ = modelEvalDescription;
}

inline
void ModelEvaluator::InArgs::_set_Np(int Np)
{
  p_.resize(Np);
}

//
// ModelEvaluator::OutArgs
//

inline
std::string ModelEvaluator::OutArgs::modelEvalDescription() const
{ return modelEvalDescription_; }

inline
int ModelEvaluator::OutArgs::Np() const
{
  return DfDp_.size();
}

inline
int ModelEvaluator::OutArgs::Ng() const
{ 
  return g_.size();
}

inline
void ModelEvaluator::OutArgs::set_f( const Evaluation<Epetra_Vector> &f ) { f_ = f; }

inline
ModelEvaluator::Evaluation<Epetra_Vector>
ModelEvaluator::OutArgs::get_f() const { return f_; }

inline
void ModelEvaluator::OutArgs::set_g( int j, const Evaluation<Epetra_Vector> &g_j )
{
  assert_j(j);
  g_[j] = g_j;
}

inline
ModelEvaluator::Evaluation<Epetra_Vector>
ModelEvaluator::OutArgs::get_g(int j) const
{
  assert_j(j);
  return g_[j];
}

inline
void ModelEvaluator::OutArgs::set_W( const Teuchos::RefCountPtr<Epetra_Operator> &W ) { W_ = W; }

inline
Teuchos::RefCountPtr<Epetra_Operator> ModelEvaluator::OutArgs::get_W() const { return W_; }

inline
ModelEvaluator::DerivativeProperties ModelEvaluator::OutArgs::get_W_properties() const
{
  return W_properties_;
}

inline
void ModelEvaluator::OutArgs::set_DfDp( int l, const Derivative &DfDp_l )
{
  assert_supports(OUT_ARG_DfDp,l);
  DfDp_[l] = DfDp_l;
}

inline
ModelEvaluator::Derivative
ModelEvaluator::OutArgs::get_DfDp(int l) const
{
  assert_supports(OUT_ARG_DfDp,l);
  return DfDp_[l];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DfDp_properties(int l) const
{
  assert_supports(OUT_ARG_DfDp,l);
  return DfDp_properties_[l];
}

inline
void ModelEvaluator::OutArgs::set_DgDx( int j, const Derivative &DgDx_j )
{
  assert_supports(OUT_ARG_DgDx,j);
  DgDx_[j] = DgDx_j;
}

inline
ModelEvaluator::Derivative
ModelEvaluator::OutArgs::get_DgDx(int j) const
{
  assert_supports(OUT_ARG_DgDx,j);
  return DgDx_[j];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DgDx_properties(int j) const
{
  assert_supports(OUT_ARG_DgDx,j);
  return DgDx_properties_[j];
}

inline
void ModelEvaluator::OutArgs::set_DgDp( int j, int l, const Derivative &DgDp_j_l )
{
  assert_supports(OUT_ARG_DgDp,j,l);
  DgDp_[ j*Np() + l ] = DgDp_j_l;
}

inline
ModelEvaluator::Derivative
ModelEvaluator::OutArgs::get_DgDp(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp,j,l);
  return DgDp_[ j*Np() + l ];
}

inline
ModelEvaluator::DerivativeProperties
ModelEvaluator::OutArgs::get_DgDp_properties(int j, int l) const
{
  assert_supports(OUT_ARG_DgDp,j,l);
  return DgDp_properties_[ j*Np() + l ];
}

inline
void ModelEvaluator::OutArgs::set_f_poly( const Teuchos::RefCountPtr<Teuchos::Polynomial<Epetra_Vector> > &f_poly )
{ f_poly_ = f_poly; }

inline
Teuchos::RefCountPtr<Teuchos::Polynomial<Epetra_Vector> >
ModelEvaluator::OutArgs::get_f_poly() const
{ return f_poly_; }

//
// ModelEvaluator::InArgsSetup
//

inline
void ModelEvaluator::InArgsSetup::setModelEvalDescription( const std::string &modelEvalDescription )
{
  this->_setModelEvalDescription(modelEvalDescription);
}

inline
void ModelEvaluator::InArgsSetup::set_Np(int Np)
{ this->_set_Np(Np); }

inline
void ModelEvaluator::InArgsSetup::setSupports( EInArgsMembers arg, bool supports )
{ this->_setSupports(arg,supports); }

//
// ModelEvaluator::OutArgsSetup
//

inline
void ModelEvaluator::OutArgsSetup::setModelEvalDescription( const std::string &modelEvalDescription )
{
  this->_setModelEvalDescription(modelEvalDescription);
}

inline
void ModelEvaluator::OutArgsSetup::set_Np_Ng(int Np, int Ng)
{ this->_set_Np_Ng(Np,Ng); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsMembers arg, bool supports )
{ this->_setSupports(arg,supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDfDp arg, int l, const DerivativeSupport& supports )
{ this->_setSupports(arg,l,supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDgDx arg, int j, const DerivativeSupport& supports )
{ this->_setSupports(arg,j,supports); }

inline
void ModelEvaluator::OutArgsSetup::setSupports( EOutArgsDgDp arg, int j, int l, const DerivativeSupport& supports )
{ this->_setSupports(arg,j,l,supports); }

inline
void ModelEvaluator::OutArgsSetup::set_W_properties( const DerivativeProperties &properties )
{ this->_set_W_properties(properties); }

inline
void ModelEvaluator::OutArgsSetup::set_DfDp_properties( int l, const DerivativeProperties &properties )
{
  this->_set_DfDp_properties(l,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DgDx_properties( int j, const DerivativeProperties &properties )
{
  this->_set_DgDx_properties(j,properties);
}

inline
void ModelEvaluator::OutArgsSetup::set_DgDp_properties( int j, int l, const DerivativeProperties &properties )
{
  this->_set_DgDp_properties(j,l,properties);
}

} // namespace EpetraExt

#endif // EPETRA_EXT_MODEL_EVALUATOR_HPP
