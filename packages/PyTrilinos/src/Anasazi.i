// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
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
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

// This documentation string will be the python help facility help
// string
%define %anasazi_docstring
"
PyTrilinos.Anasazi is the python interface to Trilinos package
Anasazi:

    http://trilinos.sandia.gov/packages/anasazi

Anasazi is a collection of eigensolver technologies.  Currently, this
python module is a stub, and none of the Anasazi classes are yet
supported.
"
%enddef

// Define the module name, its package and documentation string
%module(package   = "PyTrilinos",
	autodoc   = "1",
	docstring = %anasazi_docstring) Anasazi

%{
// System includes
#include <sstream>

// Configuration includes
#include "PyTrilinos_config.h"
#include "NumPyImporter.h"

// Teuchos includes
#include "Teuchos_ScalarTraits.hpp"

// Teuchos/python include
#include "Teuchos_PythonParameter.h"

// Epetra includes
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_SerialDenseSVD.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_SerialDistributor.h"
#include "Epetra_MapColoring.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_BasicRowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_FEVbrMatrix.h"
#include "Epetra_JadMatrix.h"

// Epetra NumPy includes
#include "Epetra_NumPyIntSerialDenseMatrix.h"
#include "Epetra_NumPyIntSerialDenseVector.h"
#include "Epetra_NumPySerialDenseMatrix.h"
#include "Epetra_NumPySerialDenseVector.h"
#include "Epetra_NumPyIntVector.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyFEVector.h"

// Anasazi includes
#include "AnasaziVersion.cpp"
#include "AnasaziTypes.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziStatusTest.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziStatusTestOrderedResNorm.hpp"
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziOrthoManager.hpp"
#include "AnasaziMatOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziEigensolverDecl.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziEpetraAdapter.hpp"

%}

// General ignore directives
%ignore *::operator=;
%ignore *::print;

// Auto-documentation feature
%feature("autodoc", "1");

// C++ STL support
%include "stl.i"

// Support for other Trilinos packages
%include "numpy.i"
%include "Teuchos_Epetra.i"

//////////////////////////////////////////////
// Support these classes, encapsulated in a //
// Teuchos::RCP<...>, as function arguments //
//////////////////////////////////////////////
%teuchos_rcp_typemaps(Epetra_MultiVector)
%teuchos_rcp_typemaps(Epetra_Operator)
%teuchos_rcp_typemaps(Anasazi::MultiVec< double >)
%teuchos_rcp_typemaps(Anasazi::StatusTest< double, Epetra_MultiVector, Epetra_Operator >)
%teuchos_rcp_typemaps(Anasazi::SortManager< double, Epetra_MultiVector, Epetra_Operator >)
%teuchos_rcp_typemaps(Anasazi::OutputManager< double >)
%teuchos_rcp_typemaps(Anasazi::Eigenproblem< double, Epetra_MultiVector, Epetra_Operator >)
%teuchos_rcp_typemaps(Anasazi::OrthoManager< double, Epetra_MultiVector >)
%teuchos_rcp_typemaps(Anasazi::MatOrthoManager< double, Epetra_MultiVector >)

/////////////////////////////////////////////////////////////////////////
// Anasazi returns eigenvalues in a std::vector< Anasazi::Value<       //
// ScalarType > > container.  This is support for converting them to a //
// more convenient numpy array.                                        //
/////////////////////////////////////////////////////////////////////////
%define %anasazi_eigenvalues_typemap(ScalarType, NumPyType)
%typemap(out)
  (std::vector< Anasazi::Value< ScalarType > >)
{
  npy_intp dims[1] = { $1.size() };
  PyObject * array = PyArray_SimpleNew(1, dims, NumPyType);
  ScalarType * data = (ScalarType*) array_data(array);
  for (npy_intp i=0; i<dims[0]; ++i) {
    data[2*i  ] = $1[i].realpart;
    data[2*i+1] = $1[i].imagpart;
  }
  return array;
  //return PyArray_SimpleNewFromData(1, dims, NumPyType, (void*)(&$1[0]));
}
%enddef
%anasazi_eigenvalues_typemap(float , NPY_CFLOAT )
%anasazi_eigenvalues_typemap(double, NPY_CDOUBLE)

/////////////////////////////
// Anasazi Version support //
/////////////////////////////
%include "AnasaziVersion.cpp"
%pythoncode %{
__version__ = Anasazi_Version().split()[2]
%}

///////////////////////////
// Anasazi Types support //
///////////////////////////
%extend Anasazi::Eigensolution {
  std::vector< Anasazi::Value< ScalarType > > Evals() {
    return self->Evals;
  }
  MV & Evecs() {
    return *(self->Evecs);
  }
  MV & Espace() {
    return *(self->Espace);
  }
}
%ignore Anasazi::Eigensolution::Evals;
%ignore Anasazi::Eigensolution::Evecs;
%ignore Anasazi::Eigensolution::Espace;
%include "AnasaziTypes.hpp"
%extend Anasazi::Value {
  std::string __str__() {
    std::stringstream output;
    output << self->realpart << "+" << self->imagpart << "j";
    return output.str();
  }
}
%template (ValueDouble)
  Anasazi::Value<double>;

///////////////////////////////////
// Anasazi OutputManager support //
///////////////////////////////////
%include "AnasaziOutputManager.hpp"
%template (OutputManagerDouble)
  Anasazi::OutputManager<double>;

////////////////////////////////////////
// Anasazi BasicOutputManager support //
////////////////////////////////////////
%include "AnasaziBasicOutputManager.hpp"
%template (BasicOutputManagerDouble)
  Anasazi::BasicOutputManager<double>;

/////////////////////////////////
// Anasazi SortManager support //
/////////////////////////////////
%include "AnasaziSortManager.hpp"

///////////////////////////////
// Anasazi BasicSort support //
///////////////////////////////
%include "AnasaziBasicSort.hpp"

////////////////////////////////////
// Anasazi MultiVecTraits support //
////////////////////////////////////
%include "AnasaziMultiVecTraits.hpp"

//////////////////////////////
// Anasazi MultiVec support //
//////////////////////////////
%include "AnasaziMultiVec.hpp"
%template (MultiVecDouble)
  Anasazi::MultiVec<double>;

////////////////////////////////////
// Anasazi OperatorTraits support //
////////////////////////////////////
%include "AnasaziOperatorTraits.hpp"

//////////////////////////////
// Anasazi Operator support //
//////////////////////////////
%include "AnasaziOperator.hpp"
%template (OperatorDouble)
  Anasazi::Operator<double>;

//////////////////////////////////
// Anasazi Eigenproblem support //
//////////////////////////////////
%include "AnasaziEigenproblem.hpp"

///////////////////////////////////////
// Anasazi BasicEigenproblem support //
///////////////////////////////////////
%include "AnasaziBasicEigenproblem.hpp"

////////////////////////////////
// Anasazi StatusTest support //
////////////////////////////////
%include "AnasaziStatusTest.hpp"

/////////////////////////////////////
// Anasazi StatusTestCombo support //
/////////////////////////////////////
%include "AnasaziStatusTestCombo.hpp"

////////////////////////////////////////
// Anasazi StatusTestMaxIters support //
////////////////////////////////////////
%include "AnasaziStatusTestMaxIters.hpp"

//////////////////////////////////////////////
// Anasazi StatusTestOrderedResNorm support //
//////////////////////////////////////////////
%include "AnasaziStatusTestOrderedResNorm.hpp"

//////////////////////////////////////
// Anasazi StatusTestOutput support //
//////////////////////////////////////
%include "AnasaziStatusTestOutput.hpp"

///////////////////////////////////////
// Anasazi StatusTestResNorm support //
///////////////////////////////////////
%include "AnasaziStatusTestResNorm.hpp"

//////////////////////////////////
// Anasazi OrthoManager support //
//////////////////////////////////
%include "AnasaziOrthoManager.hpp"

/////////////////////////////////////
// Anasazi MatOrthoManager support //
/////////////////////////////////////
%include "AnasaziMatOrthoManager.hpp"

///////////////////////////////////////
// Anasazi BasicOrthoManager support //
///////////////////////////////////////
%include "AnasaziBasicOrthoManager.hpp"

//////////////////////////////////////
// Anasazi SVQBOrthoManager support //
//////////////////////////////////////
%include "AnasaziSVQBOrthoManager.hpp"

/////////////////////////////////
// Anasazi Eigensolver support //
/////////////////////////////////
%include "AnasaziEigensolver.hpp"

///////////////////////////////////
// Anasazi SolverManager support //
///////////////////////////////////
%include "AnasaziSolverManager.hpp"

///////////////////////////////////
// Anasazi BlockDavidson support //
///////////////////////////////////
%include "AnasaziBlockDavidson.hpp"

/////////////////////////////////////////
// Anasazi BlockDavidsonSolMgr support //
/////////////////////////////////////////
%include "AnasaziBlockDavidsonSolMgr.hpp"

//////////////////////////////////////
// Anasazi BlockKrylovSchur support //
//////////////////////////////////////
%include "AnasaziBlockKrylovSchur.hpp"

////////////////////////////////////////////
// Anasazi BlockKrylovSchurSolMgr support //
////////////////////////////////////////////
%include "AnasaziBlockKrylovSchurSolMgr.hpp"

////////////////////////////
// Anasazi LOBPCG support //
////////////////////////////
%include "AnasaziLOBPCG.hpp"

//////////////////////////////////
// Anasazi LOBPCGSolMgr support //
//////////////////////////////////
%include "AnasaziLOBPCGSolMgr.hpp"

///////////////////////////////////
// Anasazi EpetraAdapter support //
///////////////////////////////////
%include "AnasaziEpetraAdapter.hpp"
%ignore
Anasazi::MultiVecTraits< double,
			 Epetra_MultiVector >::CloneView(const Epetra_MultiVector &,
							 const std::vector< int > &);
%template (SortManagerEpetra)
  Anasazi::SortManager< double, Epetra_MultiVector, Epetra_Operator >;
%template (BasicSortEpetra)
  Anasazi::BasicSort< double, Epetra_MultiVector, Epetra_Operator >;
%template (MultiVecTraitsEpetra)
  Anasazi::MultiVecTraits< double, Epetra_MultiVector >;
%template (OperatorTraitsEpetra)
  Anasazi::OperatorTraits< double, Epetra_MultiVector, Epetra_Operator >;
%template (EigenproblemEpetra)
  Anasazi::Eigenproblem< double, Epetra_MultiVector, Epetra_Operator >;
%template (BasicEigenproblemEpetra)
  Anasazi::BasicEigenproblem< double, Epetra_MultiVector, Epetra_Operator >;
%template (StatusTestEpetra)
  Anasazi::StatusTest< double, Epetra_MultiVector, Epetra_Operator >;
%template (StatusTestComboEpetra)
  Anasazi::StatusTestCombo< double, Epetra_MultiVector, Epetra_Operator >;
%template (StatusTestMaxItersEpetra)
  Anasazi::StatusTestMaxIters< double, Epetra_MultiVector, Epetra_Operator >;
%template (StatusTestOrderedResNormEpetra)
  Anasazi::StatusTestOrderedResNorm< double, Epetra_MultiVector, Epetra_Operator >;
%template (StatusTestOutputEpetra)
  Anasazi::StatusTestOutput< double, Epetra_MultiVector, Epetra_Operator >;
%template (StatusTestResNormEpetra)
  Anasazi::StatusTestResNorm< double, Epetra_MultiVector, Epetra_Operator >;
%template (OrthoManagerEpetra)
  Anasazi::OrthoManager< double, Epetra_MultiVector >;
%template (MatOrthoManagerEpetra)
  Anasazi::MatOrthoManager< double, Epetra_MultiVector, Epetra_Operator >;
%template (BasicOrthoManagerEpetra)
  Anasazi::BasicOrthoManager< double, Epetra_MultiVector, Epetra_Operator >;
%template (SVQBOrthoManagerEpetra)
  Anasazi::SVQBOrthoManager< double, Epetra_MultiVector, Epetra_Operator >;
%template (EigensolverEpetra)
  Anasazi::Eigensolver<  double, Epetra_MultiVector, Epetra_Operator >;
%template (SolverManagerEpetra)
  Anasazi::SolverManager<  double, Epetra_MultiVector, Epetra_Operator >;
%template (BlockDavidsonEpetra)
  Anasazi::BlockDavidson< double, Epetra_MultiVector, Epetra_Operator >;
%template (BlockDavidsonSolMgrEpetra)
  Anasazi::BlockDavidsonSolMgr< double, Epetra_MultiVector, Epetra_Operator >;
%template (BlockKrylovSchurEpetra)
  Anasazi::BlockKrylovSchur< double, Epetra_MultiVector, Epetra_Operator >;
%template (BlockKrylovSchurSolMgrEpetra)
  Anasazi::BlockKrylovSchurSolMgr< double, Epetra_MultiVector, Epetra_Operator >;
%template (LOBPCGEpetra)
  Anasazi::LOBPCG< double, Epetra_MultiVector, Epetra_Operator >;
%template (LOBPCGSolMgrEpetra)
  Anasazi::LOBPCGSolMgr< double, Epetra_MultiVector, Epetra_Operator >;
%template(EigensolutionEpetra)
  Anasazi::Eigensolution< double, Epetra_MultiVector >;

/////////////////////////
// std::vector support //
/////////////////////////
%template (VectorValueDouble) std::vector< Anasazi::Value< double > >;
%template (VectorInt        ) std::vector< int >;

//////////////////////////////
// Generic python interface //
//////////////////////////////
%define %anasazi_scalartype_factory(ClassName)
%pythoncode %{
def ClassName(*args):
    return ClassName##Double(*args)
%}
%enddef

%define %anasazi_factory(ClassName)
%pythoncode %{
def ClassName(*args):
    return ClassName##Epetra(*args)
%}
%enddef

// Implement the python factories for classes templated purely on ScalarType
%anasazi_scalartype_factory(Value             )
%anasazi_scalartype_factory(OutputManager     )
%anasazi_scalartype_factory(BasicOutputManager)
%anasazi_scalartype_factory(MultiVec          )
%anasazi_scalartype_factory(Operator          )

// Implement the python factories for classes templated on MV and/or OP 
%anasazi_factory(SortManager             )
%anasazi_factory(BasicSort               )
%anasazi_factory(MultiVecTraits          )
%anasazi_factory(OperatorTraits          )
%anasazi_factory(Eigenproblem            )
%anasazi_factory(BasicEigenproblem       )
%anasazi_factory(StatusTest              )
%anasazi_factory(StatusTestCombo         )
%anasazi_factory(StatusTestMaxIters      )
%anasazi_factory(StatusTestOrderedResNorm)
%anasazi_factory(StatusTestOutput        )
%anasazi_factory(StatusTestResNorm       )
%anasazi_factory(OrthoManager            )
%anasazi_factory(MatOrthoManager         )
%anasazi_factory(SVQBOrthoManager        )
%anasazi_factory(Eigensolver             )
%anasazi_factory(SolverManager           )
%anasazi_factory(BlockDavidson           )
%anasazi_factory(BlockDavidsonSolMgr     )
%anasazi_factory(BlockKrylovSchur        )
%anasazi_factory(BlockKrylovSchurSolMgr  )
%anasazi_factory(LOBPCG                  )
%anasazi_factory(LOBPCGSolMgr            )
%anasazi_factory(Eigensolution           )
