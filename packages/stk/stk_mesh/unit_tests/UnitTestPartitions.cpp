#include <stdexcept>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_util/environment/WallTime.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/fixtures/SelectorFixture.hpp>

// Borrow a lot from UnitTestSelector.  Bulk up the SelectorFixture to have parts
// with enough entities so that each partition (bucket family) comprises multiple
// buckets.

namespace {

using stk::mesh::fixtures::SelectorFixture ;


void initialize(SelectorFixture& fixture,
                std::vector<stk::mesh::Entity> &ec1_arg,
                std::vector<stk::mesh::Entity> &ec2_arg,
                std::vector<stk::mesh::Entity> &ec3_arg,
                std::vector<stk::mesh::Entity> &ec4_arg,
                std::vector<stk::mesh::Entity> &ec5_arg
                )
{
    fixture.m_meta_data.commit();
    fixture.m_bulk_data.modification_begin();
    fixture.generate_mesh();

    const size_t bucket_size = 1000;     // Default value for BucketRepository constructor.
    const size_t lb_num_buckets_per_partition = 3;
    stk::mesh::EntityRank ent_type = 0;  // rank

    const size_t bf_size = bucket_size * lb_num_buckets_per_partition;
    stk::mesh::EntityId ent_id = 1001;   // Want to keep numerical alignment.
    std::vector<stk::mesh::Part*> partMembership;

    // Note that the loop variables start at 1 because SelectorFixture::generate_mesh() has
    // already created an Entity in each partition.

    // Entities in collection 1 are contained in PartA
    partMembership.clear();
    partMembership.push_back( & fixture.m_partA );
    for (size_t i = 1; i < bf_size; ++i)
    {
        stk::mesh::Entity ent =  fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
        ec1_arg.push_back(ent);
        ++ent_id;
    }

    // Entity2 is contained in PartA and PartB
    ++ent_id;     // For numerical alignment.
    partMembership.clear();
    partMembership.push_back( & fixture.m_partA );
    partMembership.push_back( & fixture.m_partB );
    for (size_t i = 1; i < bf_size; ++i)
    {
        stk::mesh::Entity ent =  fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
        ec2_arg.push_back(ent);
        ++ent_id;
    }

    // Entity3 is contained in PartB and PartC
    ++ent_id;
    partMembership.clear();
    partMembership.push_back( & fixture.m_partB );
    partMembership.push_back( & fixture.m_partC );
    for (size_t i = 1; i < bf_size; ++i)
    {
        stk::mesh::Entity ent =  fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
        ec3_arg.push_back(ent);
        ++ent_id;
    }

    // Entity4 is contained in PartC
    ++ent_id;
    partMembership.clear();
    partMembership.push_back( & fixture.m_partC );
    for (size_t i = 1; i < bf_size; ++i)
    {
        stk::mesh::Entity ent =  fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
        ec4_arg.push_back(ent);
        ++ent_id;
    }

    // Entity5 is not contained in any Part
    ++ent_id;
    partMembership.clear();
    for (size_t i = 1; i < bf_size; ++i)
    {
        stk::mesh::Entity ent =  fixture.m_bulk_data.declare_entity(ent_type, ent_id, partMembership);
        ec5_arg.push_back(ent);
        ++ent_id;
    }

    STKUNIT_ASSERT(fixture.m_bulk_data.modification_end());
}

void initialize_reverse_ordered(SelectorFixture& fixture,
                                std::vector<stk::mesh::Entity> &ec1_arg,
                                std::vector<stk::mesh::Entity> &ec2_arg,
                                std::vector<stk::mesh::Entity> &ec3_arg,
                                std::vector<stk::mesh::Entity> &ec4_arg,
                                std::vector<stk::mesh::Entity> &ec5_arg )
{
    initialize(fixture, ec1_arg, ec2_arg, ec3_arg, ec4_arg, ec5_arg);

    stk::mesh::impl::BucketRepository &bucket_repository =
            stk::mesh::impl::Partition::getRepository(fixture.m_bulk_data);
    bucket_repository.sync_to_partitions();

    std::vector<stk::mesh::impl::Partition *> partitions = bucket_repository.get_partitions(0);
    size_t num_partitions = partitions.size();
    STKUNIT_EXPECT_EQ(num_partitions, 5u);

    for (size_t i = 0; i < num_partitions; ++i)
    {
        stk::mesh::impl::Partition &partition = *partitions[i];
        partition.reverseEntityOrderWithinBuckets();
    }
}


void initialize_data(SelectorFixture& fix)
{
    stk::mesh::impl::BucketRepository &bucket_repository =
            stk::mesh::impl::Partition::getRepository(fix.m_bulk_data);
    bucket_repository.sync_to_partitions();

    std::vector<stk::mesh::impl::Partition *> partitions = bucket_repository.get_partitions(0);
    size_t num_partitions = partitions.size();

    for (size_t i = 0; i < num_partitions; ++i)
    {
        stk::mesh::impl::Partition &partition = *partitions[i];
        for (std::vector<stk::mesh::Bucket *>::iterator bkt_i= partition.begin(); bkt_i != partition.end(); ++bkt_i)
        {
            stk::mesh::Bucket &bkt = **bkt_i;
            double *field_data = bkt.field_data(fix.m_fieldABC, bkt[0]);
            if (!field_data)
            {
                continue;
            }

            size_t bkt_size = bkt.size();
            for (size_t k = 0; k < bkt_size; ++k)
            {
                stk::mesh::Entity curr_ent = bkt[k];
                field_data[k] = curr_ent.identifier();
            }
        }
    }
}

bool check_bucket_ptrs(const stk::mesh::Bucket &bucket)
{
    if (bucket.size() == 0 )
    {
        std::cout << "Bucket has size zero!" << std::endl;
        return false;
    }

    stk::mesh::Bucket::iterator e_i = bucket.begin();
    stk::mesh::Bucket::iterator e_e = bucket.end();
    for(; e_i != e_e; ++e_i)
    {
        if (&bucket != e_i->bucket_ptr())
            return false;
    }
    return true;
}

bool check_strictly_ordered(const stk::mesh::Bucket &bucket)
{
    if (bucket.size() == 0 )
        return true;

    stk::mesh::EntityLess eless;
    stk::mesh::Bucket::iterator e_i = bucket.begin();
    stk::mesh::Bucket::iterator e_last = bucket.end() - 1;
    for(; e_i != e_last; ++e_i)
    {
        if (!eless(*e_i, *(e_i + 1)))
            return false;
    }
    return true;
}

template <typename Data_T>
bool check_nonempty_strictly_ordered(Data_T data[], size_t sz)
{
    if (sz == 0)
        return false;

    for (size_t i = 0; i < sz - 1; ++i)
    {
        if (data[i] >= data[i + 1])
        {
            std::cout << "i = " << i << ": data[i] = " << data[i]
                      << ", data[i + 1] = " << data[i + 1] << std::endl;
            return false;
        }
    }
    return true;
}

bool check_data_consistent(stk::mesh::Field<double> &data_field, const stk::mesh::Bucket &bucket)
{
    if (bucket.size() == 0 )
        return true;

    const double *field_data = bucket.field_data(data_field, bucket[0]);
    if (field_data)
    {
        size_t num_entities = bucket.size();
        for (size_t i = 0; i < num_entities; ++i)
        {
            stk::mesh::EntityId entity_id = bucket[i].identifier();
            double val = field_data[i];
            if ((val < entity_id - 0.001) || (entity_id + 0.001 < val))
            {
                std::cout << "val = " << val << "  entity_id = " << entity_id << std::endl;
                return false;
            }
        }
    }
    return true;
}


void check_test_partition_invariant(const SelectorFixture& fix, const stk::mesh::impl::Partition &partition)
{
    const std::vector<unsigned> &partition_key = partition.get_legacy_partition_id();
    for (std::vector<stk::mesh::Bucket *>::const_iterator bkt_i= partition.begin();
            bkt_i != partition.end(); ++bkt_i)
    {
        const stk::mesh::Bucket &bkt = **bkt_i;
        STKUNIT_EXPECT_EQ(&partition, bkt.getPartition() );
        STKUNIT_EXPECT_TRUE(check_bucket_ptrs(bkt));
        double *field_data = bkt.field_data(fix.m_fieldABC, bkt[0]);
        if (field_data)
        {
            STKUNIT_EXPECT_TRUE(check_nonempty_strictly_ordered(field_data, bkt.size()));
        }
        STKUNIT_EXPECT_TRUE(check_data_consistent(fix.m_fieldABC, bkt));
        const unsigned *bucket_key = bkt.key();
        for (size_t k = 0; k < partition_key.size() - 1; ++k)
        {
            STKUNIT_EXPECT_EQ(partition_key[k], bucket_key[k]);
        }
    }
}


/** \defgroup stk_mesh_partition_unit "stk::mesh::Partition Unit Testing"
  * \addtogroup stk_mesh_partition_unit
  * \{
  *
  * Selector unit testing environment. <br>
  * A special set of mesh parts and entities are set up in the
  * following configuration for the Selector unit testing.<br>
  * Parts:  PartA, PartB, PartC, PartD, PartU <br>
  * PartU = MetaData.universal_part() <br>
  * Entities:  Entity1, Entity2, Entity3, Entity4, Entity5 <br>
  *
  * PartA contains Entity1, Entity2 <br>
  * PartB contains Entity2, Entity3 <br>
  * PartC contains Entity3, Entity4 <br>
  * PartD contains no entities <br>
  * Entity5 is not contained in any Part <br>
  *
  * <PRE>
  * |----------|--|-------|--|----------|    |-------------|
  * |<--PartA---->|       |<--PartC---->|    |   PartD     |
  * |          |<---PartB--->|          |    |             |
  * |  1       |2 |       |3 |       4  | 5  |             |
  * |          |  |       |  |          |    |             |
  * |          |  |       |  |          |    |             |
  * |----------|--|-------|--|----------|    |-------------|
  * </PRE>
  *
  * Note:  The unit test names use the convention of "i" for
  * intersection, "u" for union, and "c" for complement.
  *
  * */

/** \brief Verify we can construct the unit testing fixture.
 *
 * */
STKUNIT_UNIT_TEST( UnitTestPartition, Partition_testInitialize )
{
    std::vector<stk::mesh::Entity> ec1;
    std::vector<stk::mesh::Entity> ec2;
    std::vector<stk::mesh::Entity> ec3;
    std::vector<stk::mesh::Entity> ec4;
    std::vector<stk::mesh::Entity> ec5;

    SelectorFixture fix;

    if (fix.m_bulk_data.parallel_size() > 1)
    {
        return;
    }

    initialize(fix, ec1, ec2, ec3, ec4, ec5);
    initialize_data(fix);

    stk::mesh::impl::BucketRepository &bucket_repository = stk::mesh::impl::Partition::getRepository(fix.m_bulk_data);
    bucket_repository.sync_to_partitions();

    std::vector<stk::mesh::impl::Partition *> partitions = bucket_repository.get_partitions(0);
    size_t num_partitions = partitions.size();
    STKUNIT_EXPECT_EQ(num_partitions, 5u);

    for (size_t i = 0; i < num_partitions; ++i)
    {
        const stk::mesh::impl::Partition &partition = *partitions[i];
        check_test_partition_invariant(fix, partition);
    }

    stk::mesh::Selector selector;

    {
        selector = fix.m_partA & !fix.m_partB;
        for (size_t i = 0; i < ec1.size(); ++i)
        {
            const stk::mesh::Bucket & bucket = ec1[i].bucket() ;
            bool result = selector(bucket);
            STKUNIT_EXPECT_TRUE(result);
        }
    }
    {
        selector = fix.m_partA & fix.m_partB;
        for (size_t i = 0; i < ec2.size(); ++i)
        {
            const stk::mesh::Bucket & bucket = ec2[i].bucket() ;
            bool result = selector(bucket);
            STKUNIT_EXPECT_TRUE(result);
        }
    }
    {
        selector = fix.m_partB & fix.m_partC;
        for (size_t i = 0; i < ec3.size(); ++i)
        {
            const stk::mesh::Bucket & bucket = ec3[i].bucket() ;
            bool result = selector(bucket);
          STKUNIT_EXPECT_TRUE(result);
        }
    }
    {
        selector = !fix.m_partB & fix.m_partC;
        for (size_t i = 0; i < ec4.size(); ++i)
        {
            const stk::mesh::Bucket & bucket = ec4[i].bucket() ;
            bool result = selector(bucket);
            STKUNIT_EXPECT_TRUE(result);
        }
    }
    {
      selector = !(fix.m_partA | fix.m_partB | fix.m_partC | fix.m_partD);
      for (size_t i = 0; i < ec5.size(); ++i)
      {
          const stk::mesh::Bucket & bucket = ec5[i].bucket() ;
          bool result = selector(bucket);
          STKUNIT_EXPECT_TRUE(result);
      }
    }
}


STKUNIT_UNIT_TEST( UnitTestPartition, Partition_getPartitions)
{
    std::vector<stk::mesh::Entity> ec1;
    std::vector<stk::mesh::Entity> ec2;
    std::vector<stk::mesh::Entity> ec3;
    std::vector<stk::mesh::Entity> ec4;
    std::vector<stk::mesh::Entity> ec5;

    SelectorFixture fix;

    if (fix.m_bulk_data.parallel_size() > 1)
    {
        return;
    }
    initialize(fix, ec1, ec2, ec3, ec4, ec5);
    initialize_data(fix);

    stk::mesh::impl::BucketRepository &bucket_repository = stk::mesh::impl::Partition::getRepository(fix.m_bulk_data);
    bucket_repository.sync_to_partitions();

    std::vector<stk::mesh::impl::Partition *> partitions = bucket_repository.get_partitions(0);
    size_t num_partitions = partitions.size();
    STKUNIT_EXPECT_EQ(num_partitions, 5u);

    for (size_t i = 0; i < num_partitions; ++i)
    {
        const stk::mesh::impl::Partition &partition = *partitions[i];
        check_test_partition_invariant(fix, partition);
    }
}


STKUNIT_UNIT_TEST( UnitTestPartition, Partition_testCompress)
{
    std::vector<stk::mesh::Entity> ec1;
    std::vector<stk::mesh::Entity> ec2;
    std::vector<stk::mesh::Entity> ec3;
    std::vector<stk::mesh::Entity> ec4;
    std::vector<stk::mesh::Entity> ec5;

    SelectorFixture fix;

    if (fix.m_bulk_data.parallel_size() > 1)
    {
        return;
    }
    initialize(fix, ec1, ec2, ec3, ec4, ec5);
    initialize_data(fix);

    stk::mesh::impl::BucketRepository &bucket_repository = stk::mesh::impl::Partition::getRepository(fix.m_bulk_data);
    bucket_repository.sync_to_partitions();

    std::vector<stk::mesh::impl::Partition *> partitions = bucket_repository.get_partitions(0);
    size_t num_partitions = partitions.size();
    STKUNIT_EXPECT_EQ(num_partitions, 5u);

    for (size_t i = 0; i < num_partitions; ++i)
    {
        stk::mesh::impl::Partition &partition = *partitions[i];
        partition.compress();
        check_test_partition_invariant(fix, partition);
    }
}


STKUNIT_UNIT_TEST( UnitTestPartition, Partition_testSort)
{
    std::vector<stk::mesh::Entity> ec1;
    std::vector<stk::mesh::Entity> ec2;
    std::vector<stk::mesh::Entity> ec3;
    std::vector<stk::mesh::Entity> ec4;
    std::vector<stk::mesh::Entity> ec5;

    SelectorFixture fix;

    if (fix.m_bulk_data.parallel_size() > 1)
    {
        return;
    }
    initialize_reverse_ordered(fix, ec1, ec2, ec3, ec4, ec5);
    initialize_data(fix);

    stk::mesh::impl::BucketRepository &bucket_repository = stk::mesh::impl::Partition::getRepository(fix.m_bulk_data);
    bucket_repository.sync_to_partitions();

    std::vector<stk::mesh::impl::Partition *> partitions = bucket_repository.get_partitions(0);
    size_t num_partitions = partitions.size();
    STKUNIT_EXPECT_EQ(num_partitions, 5u);

    for (size_t i = 0; i < num_partitions; ++i)
    {
        stk::mesh::impl::Partition &partition = *partitions[i];
        // check_test_partition_invariant(fix, partition);
        partition.sort();
        check_test_partition_invariant(fix, partition);
    }
}

STKUNIT_UNIT_TEST( UnitTestPartition, Partition_testRemove)
{
    std::vector<stk::mesh::Entity> ec1;
    std::vector<stk::mesh::Entity> ec2;
    std::vector<stk::mesh::Entity> ec3;
    std::vector<stk::mesh::Entity> ec4;
    std::vector<stk::mesh::Entity> ec5;

    SelectorFixture fix;

    if (fix.m_bulk_data.parallel_size() > 1)
    {
        return;
    }
    initialize(fix, ec1, ec2, ec3, ec4, ec5);
    initialize_data(fix);

    stk::mesh::impl::BucketRepository &bucket_repository = stk::mesh::impl::Partition::getRepository(fix.m_bulk_data);
    bucket_repository.sync_to_partitions();

    std::vector<stk::mesh::impl::Partition *> partitions = bucket_repository.get_partitions(0);
    size_t num_partitions = partitions.size();
    STKUNIT_EXPECT_EQ(num_partitions, 5u);

    for (size_t i = 0; i < num_partitions; ++i)
    {
        stk::mesh::impl::Partition &partition = *partitions[i];

        // Remove non-last entity in a bucket.
        stk::mesh::Bucket &bkt_0 = **partition.begin();
        partition.remove(bkt_0[0]);

        // Remove last entity in a bucket.
        stk::mesh::Bucket &bkt_1 = **partition.begin();
        stk::mesh::Entity e_last_in_1 = bkt_1[bkt_1.size() - 1];
        partition.remove(e_last_in_1);

        // Empty out the last bucket.
        stk::mesh::Bucket *last_bkt = *(partition.end() - 1);
        while (last_bkt == *(partition.end() - 1))
        {
            partition.remove((*last_bkt)[0]);
        }

        partition.sort();

        check_test_partition_invariant(fix, partition);
    }
}


/** \} */


} // namespace
