// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_TILED_CRS_PRODUCT_TENSOR_HPP
#define STOKHOS_CUDA_TILED_CRS_PRODUCT_TENSOR_HPP

#include <iostream>

#include "KokkosArray_Cuda.hpp"
#include "Cuda/KokkosArray_Cuda_Parallel.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_TiledCrsProductTensor.hpp"

#include "cuda_profiler_api.h"

namespace Stokhos {

#if 1

//----------------------------------------------------------------------------

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< TiledCrsProductTensor< TensorScalar, KokkosArray::Cuda >,
                  MatrixScalar, KokkosArray::Cuda >,
  KokkosArray::View<VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda>,
  KokkosArray::View<VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda>,
  DefaultSparseMatOps >
{
public:

  typedef KokkosArray::Cuda  device_type ;
  typedef device_type::size_type size_type ;

  typedef TiledCrsProductTensor< TensorScalar, device_type > tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type ;
  typedef KokkosArray::View< VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda >  vector_type ;

  class ProductTensorLoop {
  public:

    const matrix_type m_A;
    const vector_type m_x;
    const vector_type m_y;
    const size_type m_block_size ;

    ProductTensorLoop( const matrix_type & A,
                       const vector_type & x,
                       const vector_type & y,
                       const size_type block_size )
      : m_A( A ), m_x( x ), m_y( y ), m_block_size(block_size) {}

    __device__
    void operator()(void) const
    {
      const size_type WarpSize = KokkosArray::Impl::CudaTraits::WarpSize;

      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();
      const size_type tile_size = m_A.block.tile_size();
      //const size_type tile_dim = dim < tile_size ? dim : tile_size;
      const size_type tile_dim = tile_size;

      VectorScalar * const sh_x_k =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      VectorScalar * const sh_x_j = sh_x_k + m_block_size*tile_dim;
      VectorScalar * const sh_A_k = sh_x_j + m_block_size*tile_dim;
      VectorScalar * const sh_A_j = sh_A_k + m_block_size*tile_dim;
      VectorScalar * const sh_y   = sh_A_j + m_block_size*tile_dim;
      volatile VectorScalar * const sh_t = sh_y + tile_dim;

      const size_type nid = blockDim.x * blockDim.y;
      const size_type tid = threadIdx.x + blockDim.x * threadIdx.y;

      // blockIdx.x == row in the deterministic (finite element) system
      const size_type iBlockEntryBeg = m_A.graph.row_map[ blockIdx.x ];
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];
      size_type numBlock = (iBlockEntryEnd-iBlockEntryBeg) / m_block_size;
      const size_type remBlock = (iBlockEntryEnd-iBlockEntryBeg) % m_block_size;
      if (remBlock > 0) ++numBlock;

      // Zero y
      for (size_type i=tid; i<dim; i+=nid) {
        m_y(i,blockIdx.x) = 0.0;
      }

      // Number of Cijk tiles
      const size_type n_tile = m_A.block.num_tiles();

      // Loop over Cijk tiles
      for (size_type tile = 0; tile<n_tile; ++tile) {

        const size_type i_offset = m_A.block.offset(tile, 0);
        const size_type j_offset = m_A.block.offset(tile, 1);
        const size_type k_offset = m_A.block.offset(tile, 2);
        const size_type i_range  = m_A.block.range(tile, 0);
        const size_type j_range  = m_A.block.range(tile, 1);
        const size_type k_range  = m_A.block.range(tile, 2);
        const size_type n_row    = m_A.block.num_rows(tile);

        // Zero y
        for (size_type i=tid; i<i_range; i+=nid) {
          sh_y[i] = 0.0;
        }

        // Loop over finite element column blocks.
        size_type iBlockEntry = iBlockEntryBeg;
        for (size_type block=0; block<numBlock;
             ++block, iBlockEntry+=m_block_size) {

          const size_type block_size =
            (block == numBlock-1 && remBlock > 0) ? remBlock : m_block_size;

          // Coalesced read blocks of X and A into shared memory
          for (size_type col=0; col<block_size; ++col) {

            const size_type iBlockColumn =
              m_A.graph.entries( iBlockEntry + col );
            const VectorScalar * const x_k = & m_x( k_offset , iBlockColumn );
            const VectorScalar * const x_j = & m_x( j_offset , iBlockColumn );
            const MatrixScalar * const A_k =
              & m_A.values( k_offset , iBlockEntry + col );
            const MatrixScalar * const A_j =
              & m_A.values( j_offset , iBlockEntry + col );

            for (size_type j=tid; j<j_range; j+=nid) {
              sh_x_j[col+j*block_size] = x_j[j];
              sh_A_j[col+j*block_size] = A_j[j];
            }
            for (size_type k=tid; k<k_range; k+=nid) {
              sh_x_k[col+k*block_size] = x_k[k];
              sh_A_k[col+k*block_size] = A_k[k];
            }

          }

          __syncthreads(); // wait for X and A to be read

          // Loop over stochastic rows in this tile
          for (size_type i=threadIdx.y; i<i_range; i+=blockDim.y) {
            VectorScalar s = 0;

            // Product tensor entries which this warp will iterate:
            const size_type lBeg = m_A.block.entry_begin(tile, i);
            const size_type lEnd = m_A.block.entry_end(tile, i);

            // Loop through sparse tensor contributions with coalesced reads.
            for (size_type l=lBeg+threadIdx.x; l<lEnd; l+=blockDim.x) {

              // Read entries from the tensor
              const int j = m_A.block.coord(l,0);
              const int k = m_A.block.coord(l,1);
              const MatrixScalar v = m_A.block.value(l);

              for (size_type col=0; col<block_size; ++col) {
                const size_type jj = col+j*block_size;
                const size_type kk = col+k*block_size;
                s += v * ( sh_A_j[jj] * sh_x_k[kk] + sh_A_k[kk] * sh_x_j[jj] );
              }

            }

            // Reduction of 'y' within 'CudaTraits::WarpSize'
            sh_t[tid] = s;
            if ( threadIdx.x + 16 < WarpSize ) sh_t[tid] += sh_t[tid+16];
            if ( threadIdx.x +  8 < WarpSize ) sh_t[tid] += sh_t[tid+ 8];
            if ( threadIdx.x +  4 < WarpSize ) sh_t[tid] += sh_t[tid+ 4];
            if ( threadIdx.x +  2 < WarpSize ) sh_t[tid] += sh_t[tid+ 2];
            if ( threadIdx.x +  1 < WarpSize ) sh_t[tid] += sh_t[tid+ 1];
            if ( threadIdx.x == 0 )
              sh_y[i] += sh_t[tid];

          }

        }

        // Wait for all threads to complete the tile
        __syncthreads();

        // Store partial sum for this tile back in global memory
        for (size_type i=tid; i<i_range; i+=nid) {
          m_y( i+i_offset , blockIdx.x ) += sh_y[i];
        }
      }

    }
  };

  //------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    const size_type row_count = A.graph.row_map.dimension_0() - 1;
    const size_type tensor_dimension = A.block.dimension();
    const size_type tile_size = A.block.tile_size();
    const size_type num_tiles = A.block.num_tiles();
    //const size_type tile_dim = std::min(tensor_dimension, tile_size);
    const size_type tile_dim = tile_size;

    const size_type nWarp = 16;
    const dim3 dBlock( KokkosArray::Impl::CudaTraits::WarpSize , nWarp , 1 );
    const dim3 dGrid( row_count , 1 , 1 );

    //const size_type shcap = KokkosArray::Impl::CudaTraits::SharedMemoryCapacity;
    const size_type block_size = 5;
    const size_type shmem =
      sizeof(VectorScalar) * ((4*block_size+1)*tile_dim + dBlock.x*dBlock.y);

#if 0

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply"
              << std::endl
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
              << "  block_size(" << block_size << ")" << std::endl
              << "  shmem(" << shmem/1024 << " kB)" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dimension(" << tensor_dimension << ")" << std::endl
              << "  tile_size(" << tile_size << ")" << std::endl
              << "  num_tiles(" << num_tiles << ")" << std::endl
             ;
#endif
    //cudaProfilerStart();
    KokkosArray::Impl::cuda_parallel_launch_local_memory<<< dGrid , dBlock , shmem >>>
      ( ProductTensorLoop( A , x , y, block_size ) );
    //cudaProfilerStop();
  }
};

#elif 0

//----------------------------------------------------------------------------

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< TiledCrsProductTensor< TensorScalar, KokkosArray::Cuda >,
                  MatrixScalar, KokkosArray::Cuda >,
  KokkosArray::View<VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda>,
  KokkosArray::View<VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda>,
  DefaultSparseMatOps >
{
public:

  typedef KokkosArray::Cuda  device_type ;
  typedef device_type::size_type size_type ;

  typedef TiledCrsProductTensor< TensorScalar, device_type > tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type ;
  typedef KokkosArray::View< VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda >  vector_type ;

  class ProductTensorLoop {
  public:

    const matrix_type m_A;
    const vector_type m_x;
    const vector_type m_y;
    const size_type m_tile;

    ProductTensorLoop( const matrix_type & A,
                       const vector_type & x,
                       const vector_type & y,
                       const size_type tile )
      : m_A( A ), m_x( x ), m_y( y ), m_tile(tile) {}

    __device__
    void operator()(void) const
    {
      const size_type WarpSize = KokkosArray::Impl::CudaTraits::WarpSize;

      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();
      const size_type tile_size = m_A.block.tile_size();

      VectorScalar * const sh_x_k =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      VectorScalar * const sh_x_j = sh_x_k + tile_size;
      VectorScalar * const sh_A_k = sh_x_j + tile_size;
      VectorScalar * const sh_A_j = sh_A_k + tile_size;
      VectorScalar * const sh_y   = sh_A_j + tile_size;
      volatile VectorScalar * const sh_t = sh_y + tile_size;

      const size_type nid = blockDim.x * blockDim.y;
      const size_type tid = threadIdx.x + blockDim.x * threadIdx.y;

      // Divide warps into 4 groups for reading x_j, x_k, A_j, and A_k
      // const size_type warps_per_group = blockDim.y / 4;
      // const size_type group_lane = threadIdx.y % warps_per_group;
      // const size_type gid = threadIdx.x + blockDim.x * group_lane;
      // const size_type ngid = warps_per_group * blockDim.x;
      // const size_type group = threadIdx.y / warps_per_group;

      const size_type i_offset = m_A.block.offset(m_tile, 0);
      const size_type j_offset = m_A.block.offset(m_tile, 1);
      const size_type k_offset = m_A.block.offset(m_tile, 2);
      const size_type i_range  = m_A.block.range(m_tile, 0);
      const size_type j_range  = m_A.block.range(m_tile, 1);
      const size_type k_range  = m_A.block.range(m_tile, 2);
      const size_type n_row    = m_A.block.num_rows(m_tile);

      // blockIdx.x == row in the deterministic (finite element) system
      const size_type iBlockEntryBeg = m_A.graph.row_map[ blockIdx.x ];
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];

      // Zero y
      for (size_type i=tid; i<i_range; i+=nid) {
        sh_y[i] = 0.0;
      }

      const size_type * __restrict cijk_j = &m_A.block.coord(0,0);
      const size_type * __restrict cijk_k = &m_A.block.coord(0,1);
      const TensorScalar * __restrict cijk_v = &m_A.block.value(0);

      // Loop over finite element column blocks.
      for (size_type iBlockEntry = iBlockEntryBeg; iBlockEntry<iBlockEntryEnd;
           ++iBlockEntry) {

        __syncthreads(); // wait for threads from previous iteration to finish

        // Coalesced read blocks of X and A into shared memory
        const size_type iBlockColumn = m_A.graph.entries( iBlockEntry );
        const VectorScalar * const x_k = & m_x( k_offset , iBlockColumn );
        const VectorScalar * const x_j = & m_x( j_offset , iBlockColumn );
        const MatrixScalar * const A_k = & m_A.values( k_offset , iBlockEntry );
        const MatrixScalar * const A_j = & m_A.values( j_offset , iBlockEntry );

        for (size_type j=tid; j<j_range; j+=nid) {
          sh_x_j[j] = x_j[j];
          sh_A_j[j] = A_j[j];
        }
        for (size_type k=tid; k<k_range; k+=nid) {
          sh_x_k[k] = x_k[k];
          sh_A_k[k] = A_k[k];
        }
        // if (group == 0)
        //   for (size_type j=gid; j<j_range; j+=ngid) sh_x_j[j] = x_j[j];
        // if (group == 1)
        //   for (size_type j=gid; j<j_range; j+=ngid) sh_A_j[j] = A_j[j];
        // if (group == 2)
        //   for (size_type k=gid; k<k_range; k+=ngid) sh_x_k[k] = x_k[k];
        // if (group == 3)
        //   for (size_type k=gid; k<k_range; k+=ngid) sh_A_k[k] = A_k[k];

        __syncthreads(); // wait for X and A to be read

        // Loop over stochastic rows in this tile
        for (size_type i=threadIdx.y; i<i_range; i+=blockDim.y) {
          VectorScalar s = 0;

          // Product tensor entries which this warp will iterate:
          const size_type lBeg = m_A.block.entry_begin(m_tile, i);
          const size_type lEnd = m_A.block.entry_end(m_tile, i);

          // Loop through sparse tensor contributions with coalesced reads.
          for (size_type l=lBeg+threadIdx.x; l<lEnd; l+=blockDim.x) {

            // Read entries from the tensor
            // const int j = m_A.block.coord(l,0);
            // const int k = m_A.block.coord(l,1);
            // const MatrixScalar v = m_A.block.value(l);
            const size_type j = cijk_j[l];
            const size_type k = cijk_k[l];
            const TensorScalar v = cijk_v[l];

            s += v * ( sh_A_j[j] * sh_x_k[k] + sh_A_k[k] * sh_x_j[j] );

          }

          // Reduction of 'y' within 'CudaTraits::WarpSize'
          sh_t[tid] = s;
          if ( threadIdx.x + 16 < WarpSize ) sh_t[tid] += sh_t[tid+16];
          if ( threadIdx.x +  8 < WarpSize ) sh_t[tid] += sh_t[tid+ 8];
          if ( threadIdx.x +  4 < WarpSize ) sh_t[tid] += sh_t[tid+ 4];
          if ( threadIdx.x +  2 < WarpSize ) sh_t[tid] += sh_t[tid+ 2];
          if ( threadIdx.x +  1 < WarpSize ) sh_t[tid] += sh_t[tid+ 1];
          if ( threadIdx.x == 0 ) sh_y[i] += sh_t[tid];

        }

      }

      // Wait for all threads to complete the tile
      __syncthreads();

      // Store partial sum for this tile back in global memory
      for (size_type i=tid; i<i_range; i+=nid) {
        m_y( i+i_offset , blockIdx.x ) += sh_y[i];
      }

    }
  };

  class Zero {
  public:
    typedef typename vector_type::value_type value_type;
    typedef typename vector_type::device_type device_type;
    typedef typename device_type::size_type size_type;

    Zero( const vector_type & x ) : m_x( x ), m_d(x.dimension_0()) {}

    KOKKOSARRAY_INLINE_FUNCTION
    void operator()( const size_type j ) const {
      for (size_type i=0; i<m_d; ++i)
        m_x(i,j) = 0.0;
    }

  private:
    const vector_type m_x;
    const size_type m_d;

  };

  //------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    const size_type row_count = A.graph.row_map.dimension_0() - 1;
    const size_type tensor_dimension = A.block.dimension();
    const size_type tile_size = A.block.tile_size();
    const size_type num_tiles = A.block.num_tiles();

    const size_type nWarp = 16;
    const dim3 dBlock( KokkosArray::Impl::CudaTraits::WarpSize , nWarp , 1 );
    const dim3 dGrid( row_count , 1 , 1 );

    const size_type shmem =
      sizeof(VectorScalar) * (5*tile_size + dBlock.x*dBlock.y);

#if 1

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply"
              << std::endl
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << ")" << std::endl
              << "  shmem(" << shmem/1024 << " kB)" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dimension(" << tensor_dimension << ")" << std::endl
              << "  tile_size(" << tile_size << ")" << std::endl
              << "  num_tiles(" << num_tiles << ")" << std::endl
             ;
#endif
    //cudaProfilerStart();
    // Zero y
    KokkosArray::parallel_for( row_count , Zero(y) );

    // Loop over Cijk tiles
    for (size_type tile = 0; tile<num_tiles; ++tile) {
      KokkosArray::Impl::cuda_parallel_launch_local_memory<<< dGrid , dBlock , shmem >>>
        ( ProductTensorLoop( A , x , y, tile ) );
    }
    //cudaProfilerStop();
  }
};

#else

//----------------------------------------------------------------------------

// #define MAX_TENSOR_OFFSETS 2000
// __constant__ KokkosArray::Cuda::size_type tensor_row_begin[MAX_TENSOR_OFFSETS];

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar >
class Multiply<
  BlockCrsMatrix< TiledCrsProductTensor< TensorScalar, KokkosArray::Cuda >,
                  MatrixScalar, KokkosArray::Cuda >,
  KokkosArray::View<VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda>,
  KokkosArray::View<VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda>,
  DefaultSparseMatOps >
{
public:

  typedef KokkosArray::Cuda  device_type ;
  typedef device_type::size_type size_type ;

  typedef TiledCrsProductTensor< TensorScalar, device_type > tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type ;
  typedef KokkosArray::View< VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda >  vector_type ;

  class ProductTensorLoop {
  public:

    const matrix_type m_A;
    const vector_type m_x;
    const vector_type m_y;

    ProductTensorLoop( const matrix_type & A,
                       const vector_type & x,
                       const vector_type & y )
      : m_A( A ), m_x( x ), m_y( y ) {}

    __device__
    void operator()(void) const
    {
      const size_type WarpSize = KokkosArray::Impl::CudaTraits::WarpSize;

      // Number of bases in the stochastic system:
      const size_type dim = m_A.block.dimension();
      const size_type tile_size = m_A.block.tile_size();
      const size_type max_num_rows = m_A.block.max_num_rows();

      const size_type nid = blockDim.x * blockDim.y * blockDim.z;
      const size_type tid =
        threadIdx.x + blockDim.x * ( threadIdx.y + blockDim.y * threadIdx.z );
      const size_type thread_lane = threadIdx.x + blockDim.x * threadIdx.y;

      // blockDim.x == FEM column block-size
      VectorScalar * const sh_x_k =
        kokkos_impl_cuda_shared_memory<VectorScalar>();
      //VectorScalar * const sh_x_j = sh_x_k + blockDim.x*tile_size;
      VectorScalar * const sh_x_j = sh_x_k;
      VectorScalar * const sh_A_k = sh_x_j + blockDim.x*tile_size;
      //VectorScalar * const sh_A_j = sh_A_k + blockDim.x*tile_size;
      VectorScalar * const sh_A_j = sh_A_k;
      VectorScalar * const sh_y   = sh_A_j + blockDim.x*tile_size;
      volatile VectorScalar * const sh_t = sh_y + tile_size;
      //volatile VectorScalar * const sh_v = sh_t + nid;
      size_type * const sh_j = (size_type*) (sh_t + nid);
      size_type * const sh_k = (size_type*) (sh_j + nid);

      // blockIdx.x == row in the deterministic (finite element) system
      // These are not coalesced, but there is nothing we can do about it
      // since we only do one row per block
      const size_type iBlockEntryBeg = m_A.graph.row_map[ blockIdx.x ];
      const size_type iBlockEntryEnd = m_A.graph.row_map[ blockIdx.x + 1 ];
      size_type numBlock = (iBlockEntryEnd-iBlockEntryBeg) / blockDim.x;
      const size_type remBlock = (iBlockEntryEnd-iBlockEntryBeg) % blockDim.x;
      if (remBlock > 0) ++numBlock;

      // Zero y
      for (size_type i=tid; i<dim; i+=nid) {
        m_y(i,blockIdx.x) = 0.0;
      }

      // Number of Cijk tiles
      const size_type n_tile = m_A.block.num_tiles();

      // Loop over Cijk tiles
      for (size_type tile = 0; tile<n_tile; ++tile) {

        // These are not coalesced
        const size_type i_offset = m_A.block.offset(tile, 0);
        const size_type j_offset = m_A.block.offset(tile, 1);
        const size_type k_offset = m_A.block.offset(tile, 2);
        const size_type i_range  = m_A.block.range(tile, 0);
        const size_type j_range  = m_A.block.range(tile, 1);
        const size_type k_range  = m_A.block.range(tile, 2);
        const size_type n_row    = m_A.block.num_rows(tile);

        // Zero y
        for (size_type i=tid; i<i_range; i+=nid) {
          sh_y[i] = 0.0;
        }

        // Loop over finite element column blocks.
        // threadIdx.x == FEM column within current block
        size_type iBlockEntry = iBlockEntryBeg;
        for (size_type block=0; block<numBlock;
             ++block, iBlockEntry+=blockDim.x) {

          const size_type block_size =
            (block == numBlock-1 && remBlock > 0) ? remBlock : blockDim.x;

          // Coalesced read blocks of X and A into shared memory
          for (size_type col=0; col<block_size; ++col) {

            // This is not a coalesced read
            const size_type iBlockColumn =
              m_A.graph.entries( iBlockEntry + col );

            const VectorScalar * const x_k = & m_x( k_offset , iBlockColumn );
            const VectorScalar * const x_j = & m_x( j_offset , iBlockColumn );
            const MatrixScalar * const A_k =
              & m_A.values( k_offset , iBlockEntry + col );
            const MatrixScalar * const A_j =
              & m_A.values( j_offset , iBlockEntry + col );
            for (size_type j=tid; j<j_range; j+=nid) {
              sh_x_j[col+j*blockDim.x] = x_j[j];
              sh_A_j[col+j*blockDim.x] = A_j[j];
            }
            // for (size_type k=tid; k<k_range; k+=nid) {
            //   sh_x_k[col+k*blockDim.x] = x_k[k];
            //   sh_A_k[col+k*blockDim.x] = A_k[k];
            // }

          }

          __syncthreads(); // wait for X and A to be read

          // Loop over stochastic rows in this tile
          for (size_type i=threadIdx.z; i<i_range; i+=blockDim.z) {
            VectorScalar s = 0;

            // Product tensor entries which this warp will iterate:
            // These are not coalesced
            const size_type lBeg = m_A.block.entry_begin(tile, i);
            const size_type lEnd = m_A.block.entry_end(tile, i);
            // const size_type lBeg = tensor_row_begin[tile*(max_num_rows+1)+i];
            // const size_type lEnd = tensor_row_begin[tile*(max_num_rows+1)+i+1];

            // Loop through sparse tensor contributions with coalesced reads.
            for (size_type l=lBeg; l<lEnd; l+=WarpSize) {
              //for (size_type l=lBeg+threadIdx.y; l<lEnd; l+=blockDim.y) {

              // Read entries from the tensor
              if (l+thread_lane < lEnd) {
                sh_j[tid] = m_A.block.coord(l+thread_lane,0);
                sh_k[tid] = m_A.block.coord(l+thread_lane,1);
                sh_t[tid] = m_A.block.value(l+thread_lane);
              }

              const int mm = (l+WarpSize<lEnd) ? WarpSize : lEnd-l;
              for (size_type m=threadIdx.y; m<mm; m+=blockDim.y) {

                if (threadIdx.x < block_size) {
                  // const int j = m_A.block.coord(l+m,0);
                  // const int k = m_A.block.coord(l+m,1);
                  // const MatrixScalar v = m_A.block.value(l+m);
                  const size_type j = sh_j[m+WarpSize*threadIdx.z];
                  const size_type k = sh_k[m+WarpSize*threadIdx.z];
                  const TensorScalar v = sh_t[m+WarpSize*threadIdx.z];
                  const size_type jj = threadIdx.x + j*blockDim.x;
                  const size_type kk = threadIdx.x + k*blockDim.x;
                  s += v * ( sh_A_j[jj] * sh_x_k[kk] + sh_A_k[kk] * sh_x_j[jj] );
                }

              }

            }

            // Reduction of 'y' within 'CudaTraits::WarpSize'
            sh_t[tid] = s;
            if ( thread_lane + 16 < WarpSize ) sh_t[tid] += sh_t[tid+16];
            if ( thread_lane +  8 < WarpSize ) sh_t[tid] += sh_t[tid+ 8];
            if ( thread_lane +  4 < WarpSize ) sh_t[tid] += sh_t[tid+ 4];
            if ( thread_lane +  2 < WarpSize ) sh_t[tid] += sh_t[tid+ 2];
            if ( thread_lane +  1 < WarpSize ) sh_t[tid] += sh_t[tid+ 1];
            if ( thread_lane == 0 ) sh_y[i] += sh_t[tid];

          }

        }

        // Wait for all threads to complete the tile
        __syncthreads();

        // Store partial sum for this tile back in global memory
        for (size_type i=tid; i<i_range; i+=nid) {
          m_y( i+i_offset , blockIdx.x ) += sh_y[i];
        }
      }

    }
  };

  //------------------------------------

  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    const size_type row_count = A.graph.row_map.dimension_0() - 1;
    const size_type tensor_dimension = A.block.dimension();
    const size_type tile_size = A.block.tile_size();
    const size_type num_tiles = A.block.num_tiles();
    const size_type max_num_rows = A.block.max_num_rows();

    const size_type warp_size = KokkosArray::Impl::CudaTraits::WarpSize;
    const size_type block_size = 8;
    const size_type nWarp = 16;
    const dim3 dBlock( block_size, warp_size / block_size, nWarp );
    const dim3 dGrid( row_count , 1 , 1 );

    const size_type shmem =
      sizeof(VectorScalar) * ( (2*block_size+1)*tile_size +
                               dBlock.x*dBlock.y*dBlock.z ) +
      sizeof(size_type) * 2 * dBlock.x*dBlock.y*dBlock.z;

    const size_type cmem = sizeof(size_type)*(max_num_rows+1)*num_tiles;

#if 1

    std::cout << "Multiply< BlockCrsMatrix< CrsProductTensor ... > >::apply"
              << std::endl
              << "  grid(" << dGrid.x << "," << dGrid.y << ")" << std::endl
              << "  block(" << dBlock.x << "," << dBlock.y << "," << dBlock.z
              << ")" << std::endl
              << "  block_size(" << block_size << ")" << std::endl
              << "  shmem(" << shmem/1024 << " kB)" << std::endl
              << "  row_count(" << row_count << ")" << std::endl
              << "  tensor_dimension(" << tensor_dimension << ")" << std::endl
              << "  tile_size(" << tile_size << ")" << std::endl
              << "  num_tiles(" << num_tiles << ")" << std::endl
              << "  cmem(" << cmem/1024 << " kB)" << std::endl
             ;
#endif

    // Copy tensor row offsets to constant device memory
    // cudaMemcpyToSymbol(tensor_row_begin, A.block.row_map_ptr(), cmem, 0,
    //                    cudaMemcpyDeviceToDevice);

    //cudaProfilerStart();
    KokkosArray::Impl::cuda_parallel_launch_local_memory<<< dGrid , dBlock , shmem >>>
      ( ProductTensorLoop( A , x , y ) );
    //cudaProfilerStop();
  }
};

#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_TILED_CRS_PRODUCT_TENSOR_HPP */
