#ifndef KOKKOS_THRUSTGPUNODE_HPP_
#define KOKKOS_THRUSTGPUNODE_HPP_

#include "Kokkos_NodeAPIConfigDefs.hpp"
#include "Kokkos_CUDANodeMemoryModel.hpp"

namespace Kokkos {

/** \brief %Kokkos node interface to the Thrust library for NVIDIA CUDA-capable GPUs.
    \ingroup kokkos_node_api
    \ingroup kokkos_cuda_support
 */
class ThrustGPUNode : public CUDANodeMemoryModel {
  public:

    ThrustGPUNode(Teuchos::ParameterList &pl);

    ~ThrustGPUNode();

    //@{ Computational methods

    template <class WDP>
    static void parallel_for(int begin, int end, WDP wdp);

    template <class WDP>
    static typename WDP::ReductionType
    parallel_reduce(int begin, int end, WDP wd);

    //@} 

  private:
    int totalMem_;
};

} // namespace Kokkos

#endif
