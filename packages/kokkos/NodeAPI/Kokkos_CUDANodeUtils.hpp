#ifndef KOKKOS_CUDANODEUTILS_HPP_
#define KOKKOS_CUDANODEUTILS_HPP_

#include <cuda.h>
#include <cuda_runtime.h>

#include "Kokkos_CUDA_util_inline_runtime.h"

#include <Teuchos_ArrayRCP.hpp>

namespace Kokkos {

  class CUDANodeDeallocator {
    public:
      static void free(void *ptr);
  };

  //! \class CUDANodeCopyBackDeallocator
  /*! \brief Allocates a segment of page-locked memory associated with CUDA
      device memory; upon deallocation, performs a copy back of the allocated host
      memory before the host memory is deallocated. */
  template <class T>
  class CUDANodeCopyBackDeallocator {
    public:
      CUDANodeCopyBackDeallocator(T *devPtr, unsigned int size);

      //! Allocate the buffer, returning a Teuchos::ArrayRCP of the requested type, with the 
      Teuchos::ArrayRCP<T> alloc()const ;

      void free(void *ptr) const;
    private:
      const unsigned int size_;
      T * const devptr_;
#ifdef HAVE_KOKKOS_DEBUG
      mutable T * originalHostPtr_;
#endif
  };

  template <class T>
  CUDANodeCopyBackDeallocator<T>::CUDANodeCopyBackDeallocator(T *devPtr, unsigned int size) 
  : devptr_(devPtr), size_(size)
  { 
#ifdef HAVE_KOKKOS_DEBUG
    originalHostPtr_ = NULL;
#endif
  }

  template <class T>
  Teuchos::ArrayRCP<T>
  CUDANodeCopyBackDeallocator<T>::alloc() const {
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION( originalHostPtr_ != NULL, std::runtime_error,
        Teuchos::typeName(*this) << "::alloc(): alloc() has already been called." );
#endif
    T *hostPtr = NULL;
    // alloc page-locked ("pinned") memory on the host
    cutilSafeCallNoSync( cudaHostAlloc( (void**)&hostPtr, size_*sizeof(T), cudaHostAllocDefault) );
#ifdef HAVE_KOKKOS_DEBUG
    // save the allocated address for debug checking
    originalHostPtr_ = hostPtr; 
#endif
    // create an ARCP<T> owning this memory, with a copy of *this for the deallocator
    const bool OwnsMem = true;
    return Teuchos::arcp<T>( hostPtr, 0, size_, *this, OwnsMem );
  }

  template <class T>
  void CUDANodeCopyBackDeallocator<T>::free(void *hostPtr) const {
#ifdef HAVE_KOKKOS_DEBUG
    TEST_FOR_EXCEPTION( hostPtr != originalHostPtr_, std::logic_error,
        Teuchos::typeName(*this) << "::free(): pointer to free does not originally allocated pointer." );
    orginalHostPtr_ = NULL;
#endif
    cutilSafeCallNoSync( cudaMemcpy( devptr_, hostPtr, size_*sizeof(T), cudaMemcpyHostToDevice) );
    cutilSafeCallNoSync( cudaFreeHost( (void**)hostPtr ) );
  }

}

#endif // KOKKOS_CUDANODEUTILS_HPP_
