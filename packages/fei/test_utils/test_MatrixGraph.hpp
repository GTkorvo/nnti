/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _test_MatrixGraph_h_
#define _test_MatrixGraph_h_

#include <fei_macros.hpp>

#include <fei_defs.h>
#include <fei_SharedPtr.hpp>

#include <test_utils/tester.hpp>

#include <test_utils/testData.hpp>

#include <fei_fwd.hpp>

/** Unit-tester for fei::MatrixGraph.
*/
class test_MatrixGraph : public tester {
 public:
  test_MatrixGraph(MPI_Comm comm);
  virtual ~test_MatrixGraph();

  const char* getName()
    {
      static const char name[] = "fei::MatrixGraph";
      return((const char*)name);
    }

  int runtests();

  int test1();

  int serialtest1();

  int test2();

  int test3();

  int test4();

  int test5();

  static fei::SharedPtr<fei::MatrixGraph> create_MatrixGraph(testData* testdata,
				int localProc, int numProcs,
				bool bothFields, bool nonSymmetric,
				const char* name,
				fei::SharedPtr<fei::VectorSpace> vectorSpacePtr,
				fei::SharedPtr<fei::Factory> factory,
                                const std::string& path,
				bool turnOnDebugOutput=false);
};

#endif // _test_MatrixGraph_h_
