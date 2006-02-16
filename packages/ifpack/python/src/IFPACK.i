// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//            PyTrilinos.IFAPCK: Python Interface to IFPACK
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

%define IFPACK_DOCSTRING
"The IFPACK module allows access to The Trilinos package IFPACK.  Note
that the 'Ifpack_' prefix has been stripped from all IFPACK objects,
but that if imported with 'from PyTrilinos import IFPACK', these
objects exist in the 'IFPACK' python namespace.  Use the python help()
facility for local documentation on classes and methods, or see the
on-line documentation for more in-depth information.

The most important classes of the IFPACK module are:
- IC
- ICT
- ILU
- ILUT
- PointRelaxation

Alternatively, one can use the factory class to create a larger variety of
preconditioners; see the Doxygen documentation of Ifpack for more details.
Finally, the following functions are avaiable in the IFPACK module:
- AnalyzeMatrix()
- AnalyzeMatrixElements()
- AnalyzeVectoElements()
- PrintSparsity()
"
%enddef

%module(package="PyTrilinos", docstring=IFPACK_DOCSTRING) IFPACK

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

#include "Ifpack_ConfigDefs.h"

// Epetra includes
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_FEVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_JadOperator.h"

// Epetra python includes
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_PyOperator.h"
#include "Epetra_PyRowMatrix.h"

// Teuchos include
#include "Teuchos_ParameterList.hpp"

// IFPACK includes
#include "Ifpack.h"
#include "Ifpack_Version.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_PointRelaxation.h"

// local includes
#include "PyTeuchos_Utils.h"

%}

// Auto-documentation feature
%feature("autodoc", "1");

%ignore operator<<(ostream &, const Ifpack_Preconditioner &);// From python, use __str__
%ignore Epetra_Version();

%rename(Factory               ) Ifpack;
%rename(Version               ) Ifpack_Version;
%rename(Preconditioner        ) Ifpack_Preconditioner;
%rename(IC                    ) Ifpack_IC;
%rename(ICT                   ) Ifpack_ICT;
%rename(ILU                   ) Ifpack_ILU;
%rename(ILUT                  ) Ifpack_ILUT;
%rename(Amesos                ) Ifpack_Amesos;
%rename(PointRelaxation       ) Ifpack_PointRelaxation;
%rename(AnalyzeMatrix         ) Ifpack_Analyze;
%rename(AnalyzeMatrixElements ) Ifpack_AnalyzeMatrixElements;
%rename(AnalyzeVectorElements ) Ifpack_AnalyzeVectorElements;
%rename(PrintSparsity         ) Ifpack_PrintSparsity;

%include "Ifpack_config.h"
%include "Ifpack_ConfigDefs.h"

// typemaps
%typemap(in) (Teuchos::ParameterList& List)
{
  $1 = CreateList($input);
  if ($1 == 0)
  {
    PyErr_SetString(PyExc_ValueError, "Expecting a dictionary");
    return NULL;
  }
}

%typemap(freearg) Teuchos::ParameterList& List
{
  delete($1);
}

using namespace std;

// Epetra interface includes
%import "Epetra.i"

// IFPACK interface includes
%include "Ifpack.h"
%include "Ifpack_Version.h"
%include "Ifpack_Utils.h"
%include "Ifpack_Preconditioner.h"
%include "Ifpack_IC.h"
%include "Ifpack_ICT.h"
%include "Ifpack_ILU.h"
%include "Ifpack_ILUT.h"
%include "Ifpack_PointRelaxation.h"

%extend Ifpack_Preconditioner
{
  string __str__() {
    stringstream os;
    os << *self;
    return os.str();
  }

  void __del__()
  {
    delete self;
  }
}

%pythoncode %{
__version__ = Version().split()[2]
%}
