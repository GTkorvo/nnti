// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//           PyTrilinos.AztecOO: Python Interface to AztecOO
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

%define AZTECOO_DOCSTRING
"""
The AztecOO module allows access to The Trilinos package AztecOO.  Note
that the 'AztecOO_' prefix has been stripped from all AztecOO objects,
but that if imported with 'from PyTrilinos import AztecOO', these
objects exist in the 'AztecOO' python namespace.  Use the python help()
facility for local documentation on classes and methods, or see the
on-line documentation for more in-depth information.

The AztecOO module cannot be used without the Epetra module. You might
want to consider the IFPACK and ML modules to extend the preconditioning
capabilities of AztecOO.


*) Brief Description

AztecOO offers a suite of Krylov accelerators (like CG and GMRES), plus
several domain decomposition preconditioners with inexact local solvers. The
overlap among the subdomains can be specified by the user; each subdomain is
assigneed to a different processor.

The most important classes of the AztecOO module is:
- AztecOO


*) Example of usage

An example of usage of AztecOO is as follows. The linear system matrix is
defined as a 5-pt Laplacian on a 2D Cartesian grid. The problem has size
10000, and the corresponding linear system is solved using CG with incomplete
Cholesky preconditioner. This example can be used in serial and parallel
environments, depending on how Trilinos was configured.

from PyTrilinos import Epetra, AztecOO, Triutils
Comm = Epetra.PyComm()
Gallery = Triutils.CrsMatrixGallery(\"laplace_2d\", Comm)
Gallery.Set(\"problem_size\", 100 * 100)
Matrix = Gallery.GetMatrix()
LHS = Gallery.GetStartingSolution()
RHS = Gallery.GetRHS()
Solver = AztecOO.AztecOO(Matrix, LHS, RHS)
Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg)
Solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_dom_decomp)
Solver.SetAztecOption(AztecOO.AZ_subdomain_solve, AztecOO.AZ_icc)
Solver.SetAztecOption(AztecOO.AZ_output, 16)
Solver.Iterate(1550, 1e-5)
"""
%enddef

%module(package="PyTrilinos", docstring=AZTECOO_DOCSTRING) AztecOO

%{
// System includes
#include <iostream>
#include <sstream>
#include <vector>

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_FEVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_JadOperator.h"

// Epetra python includes
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_PyOperator.h"
#include "Epetra_PyRowMatrix.h"

// AztecOO includes
#include "AztecOO.h"
#include "AztecOO_Version.h"

%}

// Auto-documentation feature
%feature("autodoc", "1");

// Epetra interface import
using namespace std;
%import "Epetra.i"

// Amesos interface includes
%include "AztecOO.h"
%include "AztecOO_Version.h"
%include "az_aztec_defs.h"

// Extend directives
%extend AztecOO {
  double GetStatus(string what)
  {
    const double* status = self->GetAztecStatus();
    if (what == "AZ_condnum")
      return(status[AZ_condnum]);
    return(0.0);
  }
}

// Python code
%pythoncode %{
__version__ = AztecOO_Version().split()[2]
%}
