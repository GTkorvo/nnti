/*
//@HEADER
// ************************************************************************
//
//               Pliris: Parallel Dense Solver Package
//                 Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
*/

#ifdef CBLAS 
#define XCOPY zcopy
#define XSCAL zscal
#define XAXPY zaxpy
#define IXAMAX izamax
#define XASUM dzasum
#define XDOT zdotu

#else

#define XCOPY(len,v1,s1,v2,s2) zcopy_(&len,v1,&s1,v2,&s2)
#define XSCAL(len,scal,v1,s1) zscal_(&len,&scal,v1,&s1)
#define XAXPY(len,scal,v1,s1,v2,s2) zaxpy_(&len,&scal,v1,&s1,v2,&s2)
#define IXAMAX(len,v1,s1) izamax_(&len,v1,&s1)
#define XASUM(v3,len,v1,s1) dzasum_(v3,&len,v1,&s1)
#define XDOT(v3,len,v1,s1,v2,s2) zdotu_(v3,&len,v1,&s1,v2,&s2)


#endif
