/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
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
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"


/*
 * This test checks that you can use non-owning Teuchos::RCP with pointers to
 * types that are only forward declared and not defined.
 */

namespace DummyNS {class UndefinedType;}

namespace Teuchos {
TEUCHOS_TYPE_NAME_TRAITS_BUILTIN_TYPE_SPECIALIZATION(DummyNS::UndefinedType);
} // namespace Teuchos


namespace {


using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcpFromUndefRef;
using Teuchos::RCP;

using DummyNS::UndefinedType;


TEUCHOS_UNIT_TEST( RCP, ForwardDeclaredUndefined )
{
  // This test ensures that you can declare a null RCP object to an undefined
  // type without trouble.
  RCP<UndefinedType> ut_rcp;
}


TEUCHOS_UNIT_TEST( RCP, ForwardDeclaredUndefined_rcp )
{
  // This test ensures that you can set a pointer to an undefined type without
  // trouble.  Note that this has to be a non-owning RCP otherwise there will
  // be issues with the destructor call.
  UndefinedType *ut_ptr = 0;
  RCP<UndefinedType> ut_rcp =
#if defined(HAS_TEUCHOS_GET_BASE_OBJ_VOID_PTR) 
    rcpFromUndefRef(*ut_ptr)
  // In this case, you have to use rcpFromUndefRef(...) in this case instead
  // of rcpFromRef() because the latter requires the object to be defined in
  // order to call dynamic_cast<const void*>(...) in order to get the base
  // object address needed for RCPNode tracing.
#else
    rcpFromRef(*ut_ptr)
    // In this case, you can use rcpFromRef(...) because the object's baseq
    // address will not be looked up using dynamic_cast and no deallocator
    // needing to know the object's will be compiled.
#endif
    ;
}


} // namespace
