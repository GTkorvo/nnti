// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Ifpack2_Details_OverlappingRowGraph_decl.hpp"

#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION

#include "Tpetra_ETIHelperMacros.h"
#include "Tpetra_CrsGraph.hpp"
#include "Ifpack2_Details_OverlappingRowGraph_def.hpp"

#define IFPACK2_LOCAL_INSTANT(LO,GO,N) \
  template class OverlappingRowGraph<Tpetra::CrsGraph< LO, GO, N > >;     \
  template class OverlappingRowGraph<Tpetra::RowGraph< LO, GO, N > >;     \

namespace Ifpack2 {
namespace Details {

  TPETRA_ETI_MANGLING_TYPEDEFS()

// Add missing instantiations from Ifpack2
#if defined(HAVE_KOKKOSCLASSIC_KOKKOSCOMPAT) && defined(KOKKOS_HAVE_PTHREAD)
  IFPACK2_LOCAL_INSTANT(int, int, Kokkos_Compat_KokkosThreadsWrapperNode)
#endif

} // namespace Details
} // namespace Ifpack2

#endif // HAVE_IFPACK2_EXPLICIT_INSTANTIATION
