// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.  
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses, 
//         templates : new C++ techniques 
//            for scientific computing 
// 
//********************************************************
//
//  NumericalTraits class to illustrate TRAITS
//
//********************************************************
// @HEADER

#ifndef SACADO_CACHEFAD_SFADTRAITS_HPP
#define SACADO_CACHEFAD_SFADTRAITS_HPP

#include "Sacado_Traits.hpp"
#include <sstream>

// Forward declarations
namespace Sacado {
  namespace CacheFad {
    template <typename T, int Num> class SFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to SFad types
  template <typename ValueT, int Num>
  struct Promote< CacheFad::SFad<ValueT,Num>, 
		  CacheFad::SFad<ValueT,Num> > {
    typedef CacheFad::SFad<ValueT,Num> type;
  };

  //! Specialization of %Promote to SFad types
  template <typename ValueT, int Num, typename R>
  struct Promote< CacheFad::SFad<ValueT,Num>, R > {
    typedef typename ValueType< CacheFad::SFad<ValueT,Num> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef CacheFad::SFad<value_type,Num> type;
  };

  //! Specialization of %Promote to SFad types
  template <typename L, typename ValueT, int Num>
  struct Promote< L, CacheFad::SFad<ValueT, Num> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< CacheFad::SFad<ValueT,Num> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef CacheFad::SFad<value_type,Num> type;
  };

  //! Specialization of %ScalarType to SFad types
  template <typename ValueT, int Num>
  struct ScalarType< CacheFad::SFad<ValueT,Num> > {
    typedef typename CacheFad::SFad<ValueT,Num>::ScalarT type;
  };

  //! Specialization of %ValueType to SFad types
  template <typename ValueT, int Num>
  struct ValueType< CacheFad::SFad<ValueT,Num> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to SFad types
  template <typename ValueT, int Num>
  struct IsADType< CacheFad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to SFad types
  template <typename ValueT, int Num>
  struct IsScalarType< CacheFad::SFad<ValueT,Num> > {
    static const bool value = false;
  };

  //! Specialization of %Value to SFad types
  template <typename ValueT, int Num>
  struct Value< CacheFad::SFad<ValueT,Num> > {
    typedef typename ValueType< CacheFad::SFad<ValueT,Num> >::type value_type;
    static const value_type& eval(const CacheFad::SFad<ValueT,Num>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to SFad types
  template <typename ValueT, int Num>
  struct ScalarValue< CacheFad::SFad<ValueT,Num> > {
    typedef typename ValueType< CacheFad::SFad<ValueT,Num> >::type value_type;
    typedef typename ScalarType< CacheFad::SFad<ValueT,Num> >::type scalar_type;
    static const scalar_type& eval(const CacheFad::SFad<ValueT,Num>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to SFad types
  template <typename ValueT, int Num>
  struct StringName< CacheFad::SFad<ValueT,Num> > {
    static std::string eval() { 
       std::stringstream ss;
      ss << "Sacado::CacheFad::SFad< " 
	 << StringName<ValueT>::eval() << ", " << Num << " >";
      return ss.str(); 
    }
  };

  //! Specialization of %IsEqual to SFad types
  template <typename ValueT, int Num>
  struct IsEqual< CacheFad::SFad<ValueT,Num> > {
    static bool eval(const CacheFad::SFad<ValueT,Num>& x, 
		     const CacheFad::SFad<ValueT,Num>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to SFad types
  template <typename ValueT, int Num>
  struct IsStaticallySized< CacheFad::SFad<ValueT,Num> > {
    static const bool value = true;
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_Fad_ScalarTraitsImp.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to SFad types
  template <typename ValueT, int Num>
  struct PromotionTraits< Sacado::CacheFad::SFad<ValueT,Num>, 
			  Sacado::CacheFad::SFad<ValueT,Num> > {
    typedef typename Sacado::Promote< Sacado::CacheFad::SFad<ValueT,Num>,
				      Sacado::CacheFad::SFad<ValueT,Num> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to SFad types
  template <typename ValueT, int Num, typename R>
  struct PromotionTraits< Sacado::CacheFad::SFad<ValueT,Num>, R > {
    typedef typename Sacado::Promote< Sacado::CacheFad::SFad<ValueT,Num>, R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to SFad types
  template <typename L, typename ValueT, int Num>
  struct PromotionTraits< L, Sacado::CacheFad::SFad<ValueT,Num> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::CacheFad::SFad<ValueT,Num> >::type 
    promote;
  };

  //! Specializtion of %Teuchos::ScalarTraits
  template <typename ValueT, int Num>
  struct ScalarTraits< Sacado::CacheFad::SFad<ValueT,Num> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::CacheFad::SFad<ValueT,Num> >
  {};

  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename ValueT, int Num>
  struct SerializationTraits<Ordinal, Sacado::CacheFad::SFad<ValueT,Num> > :
    public Sacado::Fad::SerializationTraitsImp< Ordinal, 
						Sacado::CacheFad::SFad<ValueT,Num> > 
  {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename ValueT, int Num>
  struct ValueTypeSerializer<Ordinal, Sacado::CacheFad::SFad<ValueT,Num> > :
    public Sacado::Fad::SerializerImp< Ordinal, 
				       Sacado::CacheFad::SFad<ValueT,Num>,
				       ValueTypeSerializer<Ordinal,ValueT> > 
  {
    typedef Sacado::CacheFad::SFad<ValueT,Num> FadType;
    typedef ValueTypeSerializer<Ordinal,ValueT> ValueSerializer;
    typedef Sacado::Fad::SerializerImp< Ordinal,FadType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
			Ordinal sz = 0) :
      Base(vs, sz) {}
  };
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_CACHEFAD_SFADTRAITS_HPP
