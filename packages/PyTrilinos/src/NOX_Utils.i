// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

// NOX_Utils is accessed from several places within PyTrilinos.NOX,
// and it also has a nested class issue to work around.  So I
// concentrate all of the NOX::Utils wrapper logic here.

#if SWIG_VERSION >= 0x030000
%feature("flatnested");
#else
// Handle the NOX::Utils:Fill and Sci nested classes by defining them
// exclusively for SWIG as though they were not nested.
namespace NOX
{
class Fill {
public:
  Fill(int ntimes, char ch);
  ~Fill();
  int n;
  char c;
};
%nestedworkaround Utils::Fill;

class Sci {
public:
  Sci(double val, int precision=-1);
  ~Sci();
  double d;
  int p;
};
%nestedworkaround Utils::Sci;
}
#endif

%{
#include "NOX_Utils.H"
%}

///////////////////////
// NOX Utils support //
///////////////////////
%rename(_print) NOX::Utils::print;
%include "NOX_Utils.H"

// SWIG thinks that Fill and Sci are un-nested NOX classes, so we
// need to trick the C++ compiler into understanding these so called
// un-nested NOX types.
%{
namespace NOX
{
typedef Utils::Fill Fill;
typedef Utils::Sci  Sci ;
}
%}
