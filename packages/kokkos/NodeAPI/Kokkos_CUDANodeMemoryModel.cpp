//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Kokkos_CUDANodeMemoryModel.hpp"
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_FancyOStream.hpp>

namespace Kokkos {

  CUDANodeMemoryModel::CUDANodeMemoryModel() 
  : allocSize_(0)
  {
    clearStatistics();
  }

  void CUDANodeMemoryModel::clearStatistics() {
    numCopiesD2H_ = 0;
    numCopiesH2D_ = 0;
    numCopiesD2D_ = 0;
    bytesCopiedD2H_ = 0;
    bytesCopiedH2D_ = 0;
    bytesCopiedD2D_ = 0;
  }

  void CUDANodeMemoryModel::printStatistics(const RCP< Teuchos::FancyOStream > &os) const {
    using std::setw;
    using std::endl;
    *os << Teuchos::typeName(*this) << " memory transfer statistics" << endl
        << setw(3) << ""      << setw(4) << "" << setw(14) << "Num copies"  << setw(4) << "" << setw(14) << "Bytes copied"  << endl
        << setw(3) << "D2H"   << setw(4) << "" << setw(14) << numCopiesD2H_ << setw(4) << "" << setw(14) << bytesCopiedD2H_ << endl
        << setw(3) << "H2D"   << setw(4) << "" << setw(14) << numCopiesH2D_ << setw(4) << "" << setw(14) << bytesCopiedH2D_ << endl
        << setw(3) << "D2D"   << setw(4) << "" << setw(14) << numCopiesD2D_ << setw(4) << "" << setw(14) << bytesCopiedD2D_ << endl;
  }

}
