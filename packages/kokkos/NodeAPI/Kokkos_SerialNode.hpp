#ifndef KOKKOS_SERIALNODE_HPP_
#define KOKKOS_SERIALNODE_HPP_

#include <Kokkos_StandardNodeMemoryModel.hpp>

namespace Kokkos {

class SerialNode : public StandardMemoryModel {
  public:
    SerialNode() {}
    template <class WDP>
    void parallel_for(int beg, int end, WDP wd) {
      for (int i=beg; i != end; ++i) {
        wd.execute(i);
      }
    }

    template <class WDP>
    typename WDP::ReductionType
    parallel_reduce(int begin, int end, WDP wd) {
      typename WDP::ReductionType result = wd.identity();
      for (int i=begin; i != end; ++i) {
        result = wd.reduce( result, wd.generate(i) );
      }
      return result;
    }

};

}

#endif