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

%define %nox_solver_docstring
"
PyTrilinos.NOX.Solver is the python interface to the Solver namespace
of the Trilinos package NOX:

    http://trilinos.sandia.gov/packages/nox

The purpose of NOX.Solver is to provide a solver manager class for
NOX.  NOX.Solver provides the following user-level classes:

    * Generic  - Base class for solver managers
    * Manager  - Concrete solver manager
"
%enddef

%module(package      = "PyTrilinos.NOX",
	autodoc      = "1",
	implicitconv = "1",
	docstring    = %nox_solver_docstring) Solver

%{
// Teuchos includes
#include "Teuchos_PythonParameter.h"

// NOX includes
#include "NOX_StatusTest_Generic.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_Stagnation.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Solver_Manager.H"

// Namespace flattening
using Teuchos::RCP;
%}

// Define macro for handling exceptions thrown by NOX.Solver methods and
// constructors
%define %nox_solver_exception(className,methodName)
  %exception NOX::Solver::className::methodName {
  try {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(Teuchos::Exceptions::InvalidParameterType e) {
    PyErr_SetString(PyExc_TypeError, e.what());
    SWIG_fail;
  }
  catch(Teuchos::Exceptions::InvalidParameter e) {
    PyErr_SetString(PyExc_KeyError, e.what());
    SWIG_fail;
  }
  catch(std::runtime_error e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
  catch(std::logic_error e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
}
%enddef

// General ignore directives
%ignore operator<<(ostream &, NOX::StatusTest::StatusType );
%ignore *::print(ostream& stream, int indent = 0) const;

// Rename directives
%rename(StatusTest_Generic) NOX::StatusTest::Generic;
%rename(StatusTest_None   ) NOX::StatusTest::None;

// Trilinos imports
%import "Teuchos.i"
%import "NOX.Abstract.i"
%import "NOX.StatusTest.i"

//////////////////////////////////
// NOX::Solver::Generic support //
//////////////////////////////////
%include "NOX_Solver_Generic.H"

//////////////////////////////////
// NOX::Solver::Manager support //
//////////////////////////////////
%nox_solver_exception(Manager,Manager)
%nox_solver_exception(Manager,getList)
%ignore NOX::Solver::Manager::Manager(const RCP< NOX::Abstract::Group & >,
                                   const RCP< NOX::StatusTest::Generic & >,
                                   const RCP< Teuchos::ParameterList > &);
%ignore NOX::Solver::Manager::reset(const RCP< NOX::Abstract::Group & >,
                                 const RCP< NOX::StatusTest::Generic & >,
                                 const RCP< Teuchos::ParameterList & >);
%extend NOX::Solver::Manager {

  // The Manager constructor and reset() method both take an
  // RCP<Teuchos::ParameterList > argument.  I have typemaps for
  // RCP<...> and typemaps for Teuchos::ParameterList, but I cannot
  // get a special typemap for the two combined to work.  Therefore I
  // brute force the two functions here.

  Manager(const RCP< NOX::Abstract::Group > & grp,
	  const RCP< NOX::StatusTest::Generic > & test,
	  PyObject * dict) {

    // Initialization
    int                      res    = 0;
    Teuchos::ParameterList * params = NULL;
    void                   * argp   = NULL;
    NOX::Solver::Manager   * mgr    = NULL;
    RCP<Teuchos::ParameterList> * params_rcp = NULL;
    static swig_type_info * swig_TPL_ptr = SWIG_TypeQuery("Teuchos::ParameterList *");

    // Convert third argument to a RCP<Teuchos::ParameterList> *
    if (PyDict_Check(dict)) {
      params = Teuchos::pyDictToNewParameterList(dict);
      if (params == NULL) {
	PyErr_SetString(PyExc_ValueError,
			"Python dictionary cannot be converted to ParameterList");
	goto fail;
      }
    }
    else {
      res = SWIG_ConvertPtr(dict, &argp, swig_TPL_ptr, 0);
      if (!SWIG_IsOK(res)) {
	PyErr_SetString(PyExc_TypeError,
			"Argument 3 cannot be converted to ParameterList");
	goto fail;
      }
      params = reinterpret_cast< Teuchos::ParameterList * >(argp);
    }
    params_rcp = new RCP<Teuchos::ParameterList> (params, false);

    // Construct and return a new Manager object
    mgr = new NOX::Solver::Manager(grp, test, *params_rcp);
    delete params_rcp;
    return mgr;
  fail:
    return NULL;
  }

  bool reset(const RCP< NOX::Abstract::Group > & grp,
	     const RCP< NOX::StatusTest::Generic > & test,
	     PyObject * dict) {

    // Initialization
    int                      res    = 0;
    bool                     result = false;
    Teuchos::ParameterList * params = NULL;
    void                   * argp   = NULL;
    RCP<Teuchos::ParameterList> * params_rcp = NULL;
    static swig_type_info * swig_TPL_ptr = SWIG_TypeQuery("Teuchos::ParameterList *");

    // Convert third argument to a RCP<Teuchos::ParameterList> *
    if (PyDict_Check(dict)) {
      params = Teuchos::pyDictToNewParameterList(dict);
      if (params == NULL) {
	PyErr_SetString(PyExc_ValueError,
			"Python dictionary cannot be converted to ParameterList");
	goto fail;
      }
    }
    else {
      res = SWIG_ConvertPtr(dict, &argp, swig_TPL_ptr, 0);
      if (!SWIG_IsOK(res)) {
	PyErr_SetString(PyExc_TypeError,
			"Argument 3 cannot be converted to ParameterList");
	goto fail;
      }
      params = reinterpret_cast< Teuchos::ParameterList * >(argp);
    }
    params_rcp = new RCP<Teuchos::ParameterList> (params, false);

    // Call result() method and return the result
    result = self->reset(grp, test, *params_rcp);
    delete params_rcp;
    return result;
  fail:
    return false;
  }
}
%include "NOX_Solver_Manager.H"
