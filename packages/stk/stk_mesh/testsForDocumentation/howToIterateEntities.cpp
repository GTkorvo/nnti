#include <gtest/gtest.h>
#include <string>
#include <mpi.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_topology/topology.hpp>
#include "stkMeshTestUtils.hpp"

namespace
{
TEST(StkMeshHowTo, iterateSidesetNodesMostEfficientlyForFieldDataAccess)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    // syntax creates faces for the surface on the positive 'x-side' of the 2x2x2 cube,
    // this part is given the name 'surface_1' when it is created [create_input_mesh()]
    const std::string generatedMeshSpecification = "generated:2x2x2|sideset:X";
    stkMeshIoBroker.open_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::Field<double> &temperatureField = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "temperature");
    stk::mesh::put_field_on_all_nodes(temperatureField);
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::Part &boundaryConditionPart = *stkMeshMetaData.get_part("surface_1");
    stk::mesh::Selector boundaryNodesSelector(boundaryConditionPart);

    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();
    const stk::mesh::BucketVector &boundaryNodeBuckets = stkMeshBulkData.get_buckets(stk::topology::NODE_RANK, boundaryNodesSelector);

    double prescribedTemperatureValue = 2.0;
    std::set<stk::mesh::EntityId> boundaryNodeIds;
    for (size_t bucketIndex = 0; bucketIndex < boundaryNodeBuckets.size(); ++bucketIndex)
    {
        stk::mesh::Bucket &nodeBucket = *boundaryNodeBuckets[bucketIndex];
        double *temperatureValues = stkMeshBulkData.field_data(temperatureField, nodeBucket);
        for (size_t nodeIndex = 0; nodeIndex < nodeBucket.size(); ++nodeIndex)
        {
            stk::mesh::Entity node = nodeBucket[nodeIndex];
            boundaryNodeIds.insert(stkMeshBulkData.identifier(node));
            temperatureValues[nodeIndex] = prescribedTemperatureValue;
        }
    }

    testTemperatureFieldSetCorrectly(temperatureField, prescribedTemperatureValue, boundaryNodeIds);
}

TEST(StkMeshHowTo, iterateSidesetNodesWithFieldDataAccess)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    // syntax creates faces for the surface on the positive 'x-side' of the 2x2x2 cube,
    // this part is given the name 'surface_1' when it is created [create_input_mesh()]
    const std::string generatedMeshSpecification = "generated:2x2x2|sideset:X";
    stkMeshIoBroker.open_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::Field<double> &temperatureField = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "temperature");
    stk::mesh::put_field_on_all_nodes(temperatureField);
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::Part &boundaryConditionPart = *stkMeshMetaData.get_part("surface_1");
    stk::mesh::Selector boundaryNodesSelector(boundaryConditionPart);

    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();

    stk::mesh::EntityVector nodes;
    stkMeshBulkData.get_selected_nodes(boundaryNodesSelector, nodes);

    double prescribedTemperatureValue = 2.0;
    std::set<stk::mesh::EntityId> boundaryNodeIds;

    for (size_t nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex)
    {
        boundaryNodeIds.insert(stkMeshBulkData.identifier(nodes[nodeIndex]));
        double *temperatureValues = stkMeshBulkData.field_data(temperatureField, nodes[nodeIndex]);
        *temperatureValues = prescribedTemperatureValue;
    }

    testTemperatureFieldSetCorrectly(temperatureField, prescribedTemperatureValue, boundaryNodeIds);
}
}
