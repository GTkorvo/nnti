/*--------------------------------------------------------------------*/
/*    Copyright 2001, 2002 Sandia Corporation.                        */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Copyright 2001, 2002 Sandia Corporation, Albuquerque, NM.

#include <stk_rebalance/Partition.hpp>

using namespace stk;
namespace stk {
  using namespace rebalance;
}

Partition::Partition(ParallelMachine comm) :
  comm_(comm),
  total_number_objects_(0),
  iter_initialized_(false) 
{
}

void 
Partition::replace_mesh( const std::vector<mesh::Entity *> &mesh_objects,
                          const stk::mesh::Field<double>   * nodal_coord_ref,
                          const stk::mesh::Field<double>   * elem_weight_ref)
{
  total_number_objects_ = mesh_objects.size();
}


Partition::~Partition()
{
}

unsigned 
Partition::num_elems() const 
{ 
  return total_number_objects_ ; 
}
