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

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"


namespace {


int n = 4;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "n", &n, "Number of elements in the vectors" );
}


template<class T>
std::vector<T> generatevector(const int n_in)
{
  using Teuchos::as;
  std::vector<T> a(n_in);
  for( int i = 0; i < n_in; ++i )
    a[i] = as<T>(i); // tests non-const operator[](i)
  return a;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( vector, defaultConstruct, T )
{
  using std::vector;
  using Teuchos::as;
  vector<T> a2;
  TEST_EQUALITY_CONST( as<int>(a2.size()), 0 );
  TEST_EQUALITY_CONST( a2.empty(), true );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( vector, sizedConstruct, T )
{
  using std::vector;
  using Teuchos::as;
  using Teuchos::getConst;
  typedef typename std::vector<T>::size_type size_type;
  vector<T> a(n);
  TEST_EQUALITY_CONST( a.empty(), false );
  TEST_EQUALITY( as<int>(a.size()), n );
  TEST_COMPARE( a.max_size(), >=, as<size_type>(n) );
  TEST_COMPARE( as<int>(a.capacity()), >=, n );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( vector, operatorBracket, T )
{
  using std::vector;
  using Teuchos::as;
  out << "\nTest that a[i] == i ... ";
  vector<T> a = generatevector<T>(n);;
  bool local_success = true;
  for( int i = 0; i < n; ++i ) {
    TEST_ARRAY_ELE_EQUALITY( a, i, as<T>(i) );
  }
  if (local_success) out << "passed\n";
  else success = false;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( vector, constAt, T )
{
  using std::vector;
  using Teuchos::as;
  out << "\nTest that a.at(i) == i ...\n";
  vector<T> a = generatevector<T>(n);;
  bool local_success = true;
  for( int i = 0; i < n; ++i ) {
    TEUCHOS_TEST_EQUALITY( a.at(i), as<T>(i), out, local_success );
  }
  if (local_success) out << "passed\n";
  else success = false;
}


//
// Instantiations
//


#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, defaultConstruct, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, sizedConstruct, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, operatorBracket, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( vector, constAt, T )


UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)


} // namespace
