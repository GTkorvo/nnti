/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <unit_tests/stk_utest_macros.hpp>

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/parallel/Parallel.hpp>


#include <unit_tests/UnitTestFixture.hpp>

namespace {

  using sunit::ExampleFixture;
  using sunit::getExampleBucket;
  using sunit::getPartVector;

// (BC,2) -> B
STKUNIT_UNIT_TEST( UnitTestGetBuckets, get_involved_parts_BC_2 )
{
  ExampleFixture fix ;
  std::vector<std::string> union_names;
  union_names.push_back("PartB");
  union_names.push_back("PartC");
  stk::mesh::PartVector union_parts = getPartVector(fix,union_names);

  const stk::mesh::Bucket & bucket = getExampleBucket(fix,2);

  stk::mesh::PartVector partsInvolved;
  stk::mesh::get_involved_parts(union_parts,bucket,partsInvolved);
  STKUNIT_ASSERT_TRUE(partsInvolved.size() == 1);
  STKUNIT_EXPECT_TRUE(partsInvolved[0]->name() == "PartB");
}


// (ABC,2) -> AB
STKUNIT_UNIT_TEST( UnitTestGetBuckets, get_involved_parts_ABC_2 )
{
  ExampleFixture fix ;
  std::vector<std::string> union_names;
  union_names.push_back("PartA");
  union_names.push_back("PartB");
  union_names.push_back("PartC");
  stk::mesh::PartVector union_parts = getPartVector(fix,union_names);

  const stk::mesh::Bucket & bucket = getExampleBucket(fix,2);

  stk::mesh::PartVector partsInvolved;
  stk::mesh::get_involved_parts(union_parts,bucket,partsInvolved);
  STKUNIT_ASSERT_TRUE(partsInvolved.size() == 2);
  STKUNIT_EXPECT_TRUE(partsInvolved[0]->name() == "PartA");
  STKUNIT_EXPECT_TRUE(partsInvolved[1]->name() == "PartB");
}


// (ABC,5) -> {}
STKUNIT_UNIT_TEST( UnitTestGetBuckets, get_involved_parts_ABC_5 )
{
  ExampleFixture fix ;
  std::vector<std::string> union_names;
  union_names.push_back("PartA");
  union_names.push_back("PartB");
  union_names.push_back("PartC");
  stk::mesh::PartVector union_parts = getPartVector(fix,union_names);

  const stk::mesh::Bucket & bucket = getExampleBucket(fix,5);

  stk::mesh::PartVector partsInvolved;
  stk::mesh::get_involved_parts(union_parts,bucket,partsInvolved);
  STKUNIT_ASSERT_TRUE(partsInvolved.size() == 0);
}


// ({},1) -> {}
STKUNIT_UNIT_TEST( UnitTestGetBuckets, get_involved_parts_Zero_1 )
{
  ExampleFixture fix ;
  std::vector<std::string> union_names;
  stk::mesh::PartVector union_parts = getPartVector(fix,union_names);

  const stk::mesh::Bucket & bucket = getExampleBucket(fix,1);

  stk::mesh::PartVector partsInvolved;
  stk::mesh::get_involved_parts(union_parts,bucket,partsInvolved);
  STKUNIT_ASSERT_TRUE(partsInvolved.size() == 0);
}


} // namespace 
