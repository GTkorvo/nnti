// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
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
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_CrsMatrix.hpp"
#include "Kokkos_Version.hpp"

namespace {

  using Kokkos::DefaultNode;
  using Kokkos::CrsMatrix;

  int N = 1000;

  typedef Kokkos::DefaultNode::DefaultNodeType Node;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption("test-size",&N,"Vector length for tests.");
  }

  //
  // UNIT TESTS
  // 

  // check that default constructor zeros out, for both V and MV
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsMatrix, DefaultConstructor, Scalar, Ordinal )
  {
    CrsMatrix<Scalar,Ordinal,Node> A;
    TEST_EQUALITY_CONST(A.getNumRows(), 0);
    TEST_EQUALITY_CONST(A.getNumEntries(), 0);
  }

 #define UNIT_TEST_GROUP_ORDINAL_SCALAR( ORDINAL, SCALAR ) \
       TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( CrsMatrix, DefaultConstructor, SCALAR, ORDINAL )
 
 #define UNIT_TEST_GROUP_ORDINAL( ORDINAL ) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, int) \
          UNIT_TEST_GROUP_ORDINAL_SCALAR(ORDINAL, float)
      UNIT_TEST_GROUP_ORDINAL(int)
      typedef short int ShortInt; UNIT_TEST_GROUP_ORDINAL(ShortInt)
 
}
 
