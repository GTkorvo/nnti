#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <exception>                    // for exception
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <gtest/gtest.h>
#include <string>                       // for string
#include "Ioss_Field.h"                 // for Field, etc
#include "gtest/gtest.h"                // for AssertHelper, EXPECT_THROW
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class FieldBase; } }
TEST(StkMeshIoBroker, CheckInvalidCallOrdering)
{
    const std::string outputFilename = "invalid_checks.exo";
    MPI_Comm communicator = MPI_COMM_WORLD;

    stk::io::StkMeshIoBroker stkIo(communicator);
    const std::string exodusFileName = "generated:1x1x8";
    size_t input_index = stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    stkIo.set_active_mesh(input_index);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();
    stk::mesh::Field<double> &field0 = stkMeshMetaData.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "displacement", 1);
    stk::mesh::put_field(field0, stkMeshMetaData.universal_part());
    stkIo.populate_bulk_data();

    {
      size_t results_output_index = stkIo.create_output_mesh(outputFilename, stk::io::WRITE_RESULTS);

      stk::mesh::FieldBase *field0a = stkMeshMetaData.get_field(stk::topology::NODE_RANK, "displacement");
      stkIo.add_field(results_output_index, *field0a);
      stkIo.add_global(results_output_index, "NotTooLate", "scalar", Ioss::Field::DOUBLE);

      // Global variables defined, process_output_request does not output globals...
      EXPECT_THROW(stkIo.process_output_request(results_output_index, 0.0), std::exception);

      stkIo.begin_output_step(results_output_index, 1.0);
      stkIo.write_defined_output_fields(results_output_index);
      stkIo.end_output_step(results_output_index);

      // Try to write a global field after the step has been ended.
      EXPECT_THROW(stkIo.write_global(results_output_index, "NotTooLate", 1.0), std::exception);

      // Try to add a field after output has already been done...
      EXPECT_THROW(stkIo.add_field(results_output_index, *field0a), std::exception);

      // Try to add a global field after output has already been done...
      EXPECT_THROW(stkIo.add_global(results_output_index, "TooLate", "scalar", Ioss::Field::DOUBLE), std::exception);

      // Try to set the subset selector after output mesh has already been written.
      stk::mesh::Selector selector;
      EXPECT_THROW(stkIo.set_subset_selector(results_output_index, selector), std::exception);

      // Try to set the use_nodeset_parts_for_node_fieldssubset selector after output mesh has already been written.
      EXPECT_THROW(stkIo.use_nodeset_for_part_nodes_fields(results_output_index, true), std::exception);

      // Try to call write_defined_output_fields without beginning an output step...
      EXPECT_THROW(stkIo.write_defined_output_fields(results_output_index), std::exception);

      // Try to call end_output_step before beginning an output step...
      EXPECT_THROW(stkIo.end_output_step(results_output_index), std::exception);

      // Try to call begin_output_step() without calling end_output_step().
      stkIo.begin_output_step(results_output_index, 1.0);
      EXPECT_THROW(stkIo.begin_output_step(results_output_index, 1.0), std::exception);
      stkIo.end_output_step(results_output_index);

    }

    unlink(outputFilename.c_str());
}

