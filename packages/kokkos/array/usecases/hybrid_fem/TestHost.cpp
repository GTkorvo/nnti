
#define KOKKOS_ARRAY_BOUNDS_CHECK 1

#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>
#include <BoxMeshPartition.hpp>
#include <BoxMeshFixture.hpp>
#include <TestBoxMeshFixture.hpp>
#include <Implicit.hpp>

#include <Kokkos_Host.hpp>

#include <Kokkos_Host_macros.hpp>
#include <ParallelDataMap_macros.hpp>
#include <TestBoxMeshFixture_macros.hpp>
#include <Implicit_macros.hpp>
#include <SparseLinearSystem_macros.hpp>
#include <SparseLinearSystemFill_macros.hpp>
#include <Kokkos_Clear_macros.hpp>

#include <SparseLinearSystem_Host.hpp>

//----------------------------------------------------------------------------

void test_box_partition( bool print )
{
  const size_t ghost_layer = 1 ;
  const size_t np_max = 10000 ;

  BoxType root_box ;

  root_box[0][0] = 0 ; root_box[0][1] = 100 ;
  root_box[1][0] = 0 ; root_box[1][1] = 200 ;
  root_box[2][0] = 0 ; root_box[2][1] = 300 ;

  const size_t cell_total =
    ( root_box[0][1] - root_box[0][0] ) *
    ( root_box[1][1] - root_box[1][0] ) *
    ( root_box[2][1] - root_box[2][0] );

  for ( size_t np = 2 ; np < np_max ; np = 2 * ( np + 1 ) ) {

    std::vector<BoxType> part_boxes( np );

    box_partition_rcb( root_box , part_boxes );

    size_t cell_goal = ( cell_total + np - 1 ) / np ;
    size_t cell_max = 0 ;

    for ( size_t i = 0 ; i < np ; ++i ) {
      cell_max = std::max( cell_max , count( part_boxes[i] ) );
    }

    if ( print ) {
      std::cout << std::endl
                << "box_part( " << np 
                << " ) max( " << cell_max
                << " ) goal( " << cell_goal
                << " ) ratio( " << double(cell_max) / double(cell_goal)
                << " )" << std::endl ;
    }

    const size_t nsample = std::min(np,(size_t)4);
    const size_t stride = ( np + nsample - 1 ) / nsample ;

    for ( size_t my_part = 0 ; my_part < np ; my_part += stride ) {
      BoxType             my_use_box ;
      std::vector<size_t> my_use_id_map ;
      size_t              my_count_interior ;
      size_t              my_count_owned ;
      size_t              my_count_uses ;
      std::vector<size_t> my_recv_counts ;
      std::vector<std::vector<size_t> > my_send_map ;

      size_t count_verify = 0 ;

      box_partition_maps( root_box , part_boxes , ghost_layer , my_part ,
                          my_use_box , my_use_id_map ,
                          my_count_interior ,
                          my_count_owned ,
                          my_count_uses ,
                          my_recv_counts ,
                          my_send_map );

      count_verify = my_count_owned ;

      if ( print ) {
        std::cout << "  my_part(" << my_part << ") layout { "
                  << "P" << my_part
                  << "(" << my_count_interior
                  << "," << ( my_count_owned - my_count_interior )
                  << ")" ;
      }

      for ( size_t i = 1 ; i < np ; ++i ) {
        if ( my_recv_counts[i] ) {
          count_verify += my_recv_counts[i] ;
          const size_t ip = ( my_part + i ) % np ;

          if ( print ) {
            std::cout << " P" << ip << "(" << my_recv_counts[i] << ")" ;
          }

          // Compare recv & send lists

          BoxType             ip_use_box ;
          std::vector<size_t> ip_use_id_map ;
          size_t              ip_count_interior ;
          size_t              ip_count_owned ;
          size_t              ip_count_uses ;
          std::vector<size_t> ip_recv_counts ;
          std::vector<std::vector<size_t> > ip_send_map ;

          box_partition_maps( root_box , part_boxes , ghost_layer , ip ,
                              ip_use_box , ip_use_id_map ,
                              ip_count_interior ,
                              ip_count_owned ,
                              ip_count_uses ,
                              ip_recv_counts ,
                              ip_send_map );

          // Sent by ip, received by my_part:

          const BoxType recv_send = intersect( part_boxes[ip] , my_use_box );
          const size_t recv_send_count = count( recv_send );

          const size_t j = ( my_part + np - ip ) % np ;

          if ( recv_send_count != my_recv_counts[i] ||
               recv_send_count != ip_send_map[j].size() ) {
            throw std::runtime_error( std::string("bad recv/send map") );
          }
        }
      }
      if ( print ) { std::cout << " }" << std::endl ; }

      if ( count_verify != my_count_uses ) {
        throw std::runtime_error( std::string("bad partition map") );
      }
    }
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_host( comm::Machine machine )
{
  if ( 1 == comm::size( machine ) ) {
    test_box_partition( false );
  }

  test_box_fixture<Kokkos::Host>( machine , 100 , 200 , 300 );

  HybridFEM::Implicit::driver<double,Kokkos::Host>( "Host" , machine , 3 , 4 , 1 );
}

