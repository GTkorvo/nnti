// @HEADER
// ***********************************************************************
// 
//      Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms 
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

#ifndef RTOPPACK_RTOP_NEW_T_HELPERS_HPP
#define RTOPPACK_RTOP_NEW_T_HELPERS_HPP

#include "RTOpPack_RTOpT.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace RTOpPack {

/** \brief Simple struct for a Scalar and an Index object.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
struct ScalarIndex {
  Scalar scalar;
  Index  index;
  ScalarIndex( const Scalar &_scalar, const Index &_index ) : scalar(_scalar), index(_index) {}
};

/** \brief Simple <tt>ReductTarget</tt> subclass for simple scalar objects.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ReductTargetScalar : public ReductTarget {
public:
  ReductTargetScalar( const Scalar &scalar = Teuchos::ScalarTraits<Scalar>::zero() ) : scalar_(scalar) {}
  void set( const Scalar &scalar ) { scalar_ = scalar; }
  const Scalar& get() const { return scalar_; }
private:
  Scalar scalar_;
};

/** \brief Simple <tt>ReductTarget</tt> subclass for <tt>Scalar,Index</tt> objects.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ReductTargetScalarIndex : public ReductTarget {
public:
  ReductTargetScalarIndex(
    const Scalar &scalar = Teuchos::ScalarTraits<Scalar>::zero()
    ,const Index &index  = Teuchos::ScalarTraits<Index>::zero()
    )
    :scalarIndex_(scalar,index)
    {}
  ReductTargetScalarIndex(
    const ScalarIndex<Scalar> &scalarIndex = ScalarIndex<Scalar>(Teuchos::ScalarTraits<Scalar>::zero(),Teuchos::ScalarTraits<Index>::zero())
    )
    :scalarIndex_(scalarIndex)
    {}
  void set( const ScalarIndex<Scalar> &scalarIndex ) { scalarIndex_ = scalarIndex; }
  const ScalarIndex<Scalar>& get() const { return scalarIndex_; }
private:
  ScalarIndex<Scalar> scalarIndex_;
};

/** \brief Simple base class for all reduction operators that return a simple
 * scalar reduction object.
 *
 * Subclasses have to minimally define <tt>apply_op()</tt>.
 * Subclasses should also override <tt>reduce_reduct_objs()</tt> if
 * the reduction is not a simple summation.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ROpScalarReductionBase : virtual public RTOpT<Scalar> {
public:
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  /** \brief . */
  ROpScalarReductionBase( const Scalar &initReductObjValue = Teuchos::ScalarTraits<Scalar>::zero() )
    :RTOpT<Scalar>(""), initReductObjValue_(initReductObjValue) 
    {}
  /** \brief . */
  const Scalar& getRawVal( const ReductTarget &reduct_obj ) const
    {
      using Teuchos::dyn_cast;
      return dyn_cast<const ReductTargetScalar<Scalar> >(reduct_obj).get();
    }
  /** \brief . */
  void setRawVal( const Scalar &rawVal, ReductTarget *reduct_obj ) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION( reduct_obj==NULL, std::invalid_argument, "Error!" );
#endif
      using Teuchos::dyn_cast;
      dyn_cast<ReductTargetScalar<Scalar> >(*reduct_obj).set(rawVal);
    }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void get_reduct_type_num_entries(
    int*   num_values
    ,int*  num_indexes
    ,int*  num_chars
    ) const
    {
      *num_values = Teuchos::PrimitiveTypeTraits<Scalar>::numPrimitiveObjs();
      *num_indexes = 0;
      *num_chars = 0;
    }
  /** \brief . */
  Teuchos::RefCountPtr<ReductTarget> reduct_obj_create() const
    {
      return Teuchos::rcp(new ReductTargetScalar<Scalar>(initReductObjValue()));
    }
  /** \brief Default implementation here is for a sum. */
  void reduce_reduct_objs(
    const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj
    ) const
    {
      const Scalar in_val    = getRawVal(in_reduct_obj);
      const Scalar inout_val = getRawVal(*inout_reduct_obj);
      setRawVal( in_val + inout_val, inout_reduct_obj );
    }
  /** \brief . */
  void reduct_obj_reinit( ReductTarget* reduct_obj ) const
    {
      setRawVal( initReductObjValue(), reduct_obj );
    }
  /** \brief . */
  void extract_reduct_obj_state(
    const ReductTarget        &reduct_obj
    ,int                      num_values
    ,primitive_value_type     value_data[]
    ,int                      num_indexes
    ,index_type               index_data[]
    ,int                      num_chars
    ,char_type                char_data[]
    ) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION(
        num_values==0 || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Teuchos::PrimitiveTypeTraits<Scalar>::extractPrimitiveObjs( getRawVal(reduct_obj), num_values, value_data );
    }
  /** \brief . */
  void load_reduct_obj_state(
    int                            num_values
    ,const primitive_value_type    value_data[]
    ,int                           num_indexes
    ,const index_type              index_data[]
    ,int                           num_chars
    ,const char_type               char_data[]
    ,ReductTarget                  *reduct_obj
    ) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION(
        num_values==0 || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Scalar val = Teuchos::ScalarTraits<Scalar>::zero();
      Teuchos::PrimitiveTypeTraits<Scalar>::loadPrimitiveObjs( num_values, value_data, &val );
      setRawVal( val, reduct_obj );
    }
  //@}
protected:
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, initReductObjValue )
}; // class ROpScalarReductionBase

/** \brief Base class for all reduction operators that return a
 * <tt>ScalarIndex</tt> reduction object.
 *
 * Subclasses have to minimally define <tt>apply_op()</tt> and
 * <tt>reduce_reduct_objs()</tt>.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ROpScalarIndexReductionBase : virtual public RTOpT<Scalar> {
public:
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  /** \brief . */
  ROpScalarIndexReductionBase(
    const Scalar &initScalarReductObjValue = Teuchos::ScalarTraits<Scalar>::zero()
    ,const Index  &initIndexReductObjValue = Teuchos::ScalarTraits<Index>::zero()
    )
    :RTOpT<Scalar>("")
    ,initScalarReductObjValue_(initScalarReductObjValue)
    ,initIndexReductObjValue_(initIndexReductObjValue)
    {}
  /** \brief . */
  const ScalarIndex<Scalar>& getRawVal( const ReductTarget &reduct_obj ) const
    {
      using Teuchos::dyn_cast;
      return dyn_cast<const ReductTargetScalarIndex<Scalar> >(reduct_obj).get();
    }
  /** \brief . */
  void setRawVal( const ScalarIndex<Scalar> &rawVal, ReductTarget *reduct_obj ) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION( reduct_obj==NULL, std::invalid_argument, "Error!" );
#endif
      using Teuchos::dyn_cast;
      dyn_cast<ReductTargetScalarIndex<Scalar> >(*reduct_obj).set(rawVal);
    }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void get_reduct_type_num_entries(
    int*   num_values
    ,int*  num_indexes
    ,int*  num_chars
    ) const
    {
      *num_values = Teuchos::PrimitiveTypeTraits<Scalar>::numPrimitiveObjs();
      *num_indexes = 1;
      *num_chars = 0;
    }
  /** \brief . */
  Teuchos::RefCountPtr<ReductTarget> reduct_obj_create() const
    {
      return Teuchos::rcp(
        new ReductTargetScalarIndex<Scalar>(
          ScalarIndex<Scalar>(initScalarReductObjValue(),initIndexReductObjValue())
          )
        );
    }
  /** \brief . */
  void reduct_obj_reinit( ReductTarget* reduct_obj ) const
    {
      setRawVal( ScalarIndex<Scalar>(initScalarReductObjValue(),initIndexReductObjValue()), reduct_obj );
    }
  /** \brief . */
  void extract_reduct_obj_state(
    const ReductTarget        &reduct_obj
    ,int                      num_values
    ,primitive_value_type     value_data[]
    ,int                      num_indexes
    ,index_type               index_data[]
    ,int                      num_chars
    ,char_type                char_data[]
    ) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION(
        num_values==0 || value_data==NULL || num_indexes==0 || index_data==NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      const ScalarIndex<Scalar> &scalarIndex = getRawVal(reduct_obj);
      Teuchos::PrimitiveTypeTraits<Scalar>::extractPrimitiveObjs( scalarIndex.scalar, num_values, value_data );
      index_data[0] = scalarIndex.index;
    }
  /** \brief . */
  void load_reduct_obj_state(
    int                            num_values
    ,const primitive_value_type    value_data[]
    ,int                           num_indexes
    ,const index_type              index_data[]
    ,int                           num_chars
    ,const char_type               char_data[]
    ,ReductTarget                  *reduct_obj
    ) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION(
        num_values==0 || value_data==NULL || num_indexes==0 || index_data==NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Scalar val = Teuchos::ScalarTraits<Scalar>::zero();
      Teuchos::PrimitiveTypeTraits<Scalar>::loadPrimitiveObjs( num_values, value_data, &val );
      setRawVal( ScalarIndex<Scalar>(val,index_data[0]), reduct_obj );
    }
  //@}
protected:
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, initScalarReductObjValue )
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Index, initIndexReductObjValue )
}; // class ROpScalarIndexReductionBase

/** \breif Simple base class for all reduction operators that return a simple
 * index reduction object.
 *
 * Subclasses have to minimally define <tt>apply_op()</tt>.
 * Subclasses should also override <tt>reduce_reduct_objs()</tt> if
 * the reduction is not a simple summation.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ROpIndexReductionBase : virtual public RTOpT<Scalar> {
public:
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  /** \brief . */
  ROpIndexReductionBase( const index_type &initReductObjValue = Teuchos::ScalarTraits<index_type>::zero() )
    :RTOpT<Scalar>(""), initReductObjValue_(initReductObjValue) 
    {}
  /** \brief . */
  index_type getRawVal( const ReductTarget &reduct_obj ) const
    {
      using Teuchos::dyn_cast;
      return dyn_cast<const ReductTargetScalar<index_type> >(reduct_obj).get();
    }
  /** \brief . */
  void setRawVal( const index_type &rawVal, ReductTarget *reduct_obj ) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION( reduct_obj==NULL, std::invalid_argument, "Error!" );
#endif
      using Teuchos::dyn_cast;
      dyn_cast<ReductTargetScalar<index_type> >(*reduct_obj).set(rawVal);
    }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void get_reduct_type_num_entries(
    int*   num_values
    ,int*  num_indexes
    ,int*  num_chars
    ) const
    {
      *num_values = 0;
      *num_indexes = 1;
      *num_chars = 0;
    }
  /** \brief . */
  Teuchos::RefCountPtr<ReductTarget> reduct_obj_create() const
    {
      return Teuchos::rcp(new ReductTargetScalar<index_type>(initReductObjValue()));
    }
  /// Default implementation here is for a sum
  void reduce_reduct_objs(
    const ReductTarget& _in_reduct_obj, ReductTarget* _inout_reduct_obj
    ) const
    {
      using Teuchos::dyn_cast;
      const ReductTargetScalar<index_type> &in_reduct_obj    = dyn_cast<const ReductTargetScalar<index_type> >(_in_reduct_obj); 
      ReductTargetScalar<index_type>       &inout_reduct_obj = dyn_cast<ReductTargetScalar<index_type> >(*_inout_reduct_obj); 
      inout_reduct_obj.set( inout_reduct_obj.get() + in_reduct_obj.get() );
    }
  /** \brief . */
  void reduct_obj_reinit( ReductTarget* reduct_obj ) const
    {
      setRawVal( initReductObjValue(), reduct_obj );
    }
  /** \brief . */
  void extract_reduct_obj_state(
    const ReductTarget        &reduct_obj
    ,int                      num_values
    ,primitive_value_type     value_data[]
    ,int                      num_indexes
    ,index_type               index_data[]
    ,int                      num_chars
    ,char_type                char_data[]
    ) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION(
        num_values!=0 || value_data!=NULL || num_indexes!=1 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      index_data[0] = getRawVal(reduct_obj);
    }
  /** \brief . */
  void load_reduct_obj_state(
    int                            num_values
    ,const primitive_value_type    value_data[]
    ,int                           num_indexes
    ,const index_type              index_data[]
    ,int                           num_chars
    ,const char_type               char_data[]
    ,ReductTarget                  *reduct_obj
    ) const
    {
      using Teuchos::dyn_cast;
#ifdef _DEBUG
      TEST_FOR_EXCEPTION(
        num_values!=0 || value_data!=NULL || num_indexes!=1 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      setRawVal( index_data[0], reduct_obj );
    }
  //@}
protected:
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( index_type, initReductObjValue )
}; // class ROpIndexReductionBase

/** \breif Simple base class for all transformation operators that
 * use a single piece of Scalar data.
 *
 * Subclasses have to minimally define <tt>apply_op()</tt>.
 * Subclasses should also define access functions for changing the
 * scalar data.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ROpScalarTransformationBase : virtual public RTOpT<Scalar> {
public:
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  /** \brief . */
  ROpScalarTransformationBase( const Scalar &scalarData = Teuchos::ScalarTraits<Scalar>::zero() )
    :RTOpT<Scalar>(""), scalarData_(scalarData)
    {}
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void get_op_type_num_entries(
    int*  num_values
    ,int* num_indexes
    ,int* num_chars
    ) const
    {
      *num_values = Teuchos::PrimitiveTypeTraits<Scalar>::numPrimitiveObjs();
      *num_indexes = 0;
      *num_chars = 0;
    }
  /** \brief . */
  void extract_op_state(
    int                             num_values
    ,primitive_value_type           value_data[]
    ,int                            num_indexes
    ,index_type                     index_data[]
    ,int                            num_chars
    ,char_type                      char_data[]
    ) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION(
        num_values==0 || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Teuchos::PrimitiveTypeTraits<Scalar>::extractPrimitiveObjs( scalarData_, num_values, value_data );
    }
  /** \brief . */
  void load_op_state(
    int                           num_values
    ,const primitive_value_type   value_data[]
    ,int                          num_indexes
    ,const index_type             index_data[]
    ,int                          num_chars
    ,const char_type              char_data[]
    )
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION(
        num_values==0 || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Teuchos::PrimitiveTypeTraits<Scalar>::loadPrimitiveObjs( num_values, value_data, &scalarData_ );
    }
  //@}
protected:
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, scalarData )
}; // class ROpScalarTransformationBase

/** \breif Simple base class for all transformation operators that
 * use a pair of Scalar data members.
 *
 * Subclasses have to minimally define <tt>apply_op()</tt>.
 * Subclasses should also define access functions for changing the
 * scalar data.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class ROpScalarScalarTransformationBase : virtual public RTOpT<Scalar> {
public:
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  /** \brief . */
  ROpScalarScalarTransformationBase(
    const Scalar &scalarData1 = Teuchos::ScalarTraits<Scalar>::zero()
    ,const Scalar &scalarData2 = Teuchos::ScalarTraits<Scalar>::zero()
    )
    :RTOpT<Scalar>(""), scalarData1_(scalarData1), scalarData2_(scalarData2)
    {}
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void get_op_type_num_entries(
    int*  num_values
    ,int* num_indexes
    ,int* num_chars
    ) const
    {
      *num_values = 2 * Teuchos::PrimitiveTypeTraits<Scalar>::numPrimitiveObjs();
      *num_indexes = 0;
      *num_chars = 0;
    }
  /** \brief . */
  void extract_op_state(
    int                             num_values
    ,primitive_value_type           value_data[]
    ,int                            num_indexes
    ,index_type                     index_data[]
    ,int                            num_chars
    ,char_type                      char_data[]
    ) const
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION(
        num_values==0 || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Teuchos::PrimitiveTypeTraits<Scalar>::extractPrimitiveObjs( scalarData1_, num_values/2, value_data );
      Teuchos::PrimitiveTypeTraits<Scalar>::extractPrimitiveObjs( scalarData2_, num_values/2, value_data + num_values/2 );
    }
  /** \brief . */
  void load_op_state(
    int                           num_values
    ,const primitive_value_type   value_data[]
    ,int                          num_indexes
    ,const index_type             index_data[]
    ,int                          num_chars
    ,const char_type              char_data[]
    )
    {
#ifdef _DEBUG
      TEST_FOR_EXCEPTION(
        num_values==0 || value_data==NULL || num_indexes!=0 || index_data!=NULL || num_chars!=0 || char_data!=NULL
        ,std::invalid_argument, "Error!"
        );
#endif
      Teuchos::PrimitiveTypeTraits<Scalar>::loadPrimitiveObjs( num_values/2, value_data, &scalarData1_ );
      Teuchos::PrimitiveTypeTraits<Scalar>::loadPrimitiveObjs( num_values/2, value_data+num_values/2, &scalarData2_ );
    }
  //@}
protected:
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, scalarData1 )
  /** \brief . */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, scalarData2 )
}; // class ROpScalarTransformationBase

} // namespace RTOpPack

/** \breif Use within an apply_op(...) function implemention where num_vecs==1, num_targ_vecs==0.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_1_0( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=1 || (SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==1, sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=0 || (TARG_SUB_VECS)!=NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
    ); \
  const RTOpPack::index_type   subDim  = (SUB_VECS)[0].subDim(); \
  const Scalar                 *v0_val = (SUB_VECS)[0].values(); \
  const ptrdiff_t              v0_s    = (SUB_VECS)[0].stride()

/** \breif Use within an apply_op(...) function implemention where num_vecs==2, num_targ_vecs==0.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_2_0( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=2 || (SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==1, sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=0 || (TARG_SUB_VECS)!=NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (SUB_VECS)[0].subDim() != (SUB_VECS)[1].subDim() || \
    (SUB_VECS)[0].globalOffset() != (SUB_VECS)[1].globalOffset() \
    ,IncompatibleVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
    ); \
  const RTOpPack::index_type   subDim  = (SUB_VECS)[0].subDim(); \
  const Scalar                 *v0_val = (SUB_VECS)[0].values(); \
  const ptrdiff_t              v0_s    = (SUB_VECS)[0].stride(); \
  const Scalar                 *v1_val = (SUB_VECS)[1].values(); \
  const ptrdiff_t              v1_s    = (SUB_VECS)[1].stride()

/** \breif Use within an apply_op(...) function implemention where num_vecs==0, num_targ_vecs==1.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_0_1( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=0 || (SUB_VECS)!=NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==0, sub_vecs==NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=1 || (TARG_SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==1, targ_sub_vecs!=NULL" \
    ); \
  const RTOpPack::index_type   subDim  = (TARG_SUB_VECS)[0].subDim(); \
  Scalar                       *z0_val = (TARG_SUB_VECS)[0].values(); \
  const ptrdiff_t              z0_s    = (TARG_SUB_VECS)[0].stride()

/** \breif Use within an apply_op(...) function implemention where num_vecs==1, num_targ_vecs==1.
 */
#define RTOP_APPLY_OP_1_1( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=1 || (SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==1, sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=1 || (TARG_SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==1, targ_sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (SUB_VECS)[0].subDim() != (TARG_SUB_VECS)[0].subDim() || \
    (SUB_VECS)[0].globalOffset() != (TARG_SUB_VECS)[0].globalOffset() \
    ,IncompatibleVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
    ); \
  const RTOpPack::index_type   subDim  = (SUB_VECS)[0].subDim(); \
  const Scalar                 *v0_val = (SUB_VECS)[0].values(); \
  const ptrdiff_t              v0_s    = (SUB_VECS)[0].stride(); \
  Scalar                       *z0_val = (TARG_SUB_VECS)[0].values(); \
  const ptrdiff_t              z0_s    = (TARG_SUB_VECS)[0].stride()

/** \breif Use within an apply_op(...) function implemention where num_vecs==2, num_targ_vecs==1.
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
#define RTOP_APPLY_OP_2_1( NUM_VECS, SUB_VECS, NUM_TARG_VECS, TARG_SUB_VECS ) \
  TEST_FOR_EXCEPTION( \
    (NUM_VECS)!=2 || (SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumVecs \
    ,"Error, num_vecs="<<(NUM_VECS)<<" not allowed, only num_vecs==2, sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (NUM_TARG_VECS)!=1 || (TARG_SUB_VECS)==NULL \
    ,RTOpPack::InvalidNumTargVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==1, targ_sub_vecs!=NULL" \
    ); \
  TEST_FOR_EXCEPTION( \
    (SUB_VECS)[0].subDim() != (SUB_VECS)[1].subDim() || \
    (SUB_VECS)[0].subDim() != (TARG_SUB_VECS)[0].subDim() || \
    (SUB_VECS)[0].globalOffset() != (SUB_VECS)[1].globalOffset() || \
    (SUB_VECS)[0].globalOffset() != (TARG_SUB_VECS)[0].globalOffset() \
    ,IncompatibleVecs \
    ,"Error, num_targ_vecs="<<(NUM_TARG_VECS)<<" not allowed, only num_targ_vecs==0, targ_sub_vecs==NULL" \
    ); \
  const RTOpPack::index_type   subDim  = (SUB_VECS)[0].subDim(); \
  const Scalar                 *v0_val = (SUB_VECS)[0].values(); \
  const ptrdiff_t              v0_s    = (SUB_VECS)[0].stride(); \
  const Scalar                 *v1_val = (SUB_VECS)[1].values(); \
  const ptrdiff_t              v1_s    = (SUB_VECS)[1].stride(); \
  Scalar                       *z0_val = (TARG_SUB_VECS)[0].values(); \
  const ptrdiff_t              z0_s    = (TARG_SUB_VECS)[0].stride()

#endif // RTOPPACK_RTOP_NEW_T_HELPERS_HPP
