/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_mpi.h>

#include <test_utils/LibraryFactory.hpp>

#include <fei_LibraryWrapper.hpp>

#include <snl_fei_Factory.hpp>

#include <fei_Factory_Trilinos.hpp>
#ifdef HAVE_AZTECOO
#include <fei_Aztec_LinSysCore.hpp>
#endif

#ifdef FEI_HAVE_FETI
#include <FETI_DP_FiniteElementData.h>
#endif

#ifdef FEI_HAVE_PETSC
#include <fei_PETSc_LinSysCore.hpp>
#endif

#ifdef FEI_HAVE_PROMETHEUS
#include <Prometheus_LinSysCore.h>
#endif

#include <test_utils/LinSysCore.hpp>

#ifdef FEI_HAVE_HYPRE
#include <HYPRE.hpp>
#include <HYPRE_config.hpp>
#include <IJ_mv/HYPRE_IJ_mv.hpp>
#include <parcsr_mv/HYPRE_parcsr_mv.hpp>
#include <parcsr_ls/HYPRE_parcsr_ls.hpp>
#include <HYPRE_parcsr_bicgstabl.hpp>
#include <HYPRE_parcsr_TFQmr.hpp>
#include <HYPRE_parcsr_bicgs.hpp>
#include <HYPRE_LinSysCore.hpp>
#endif

//----------------------------------------------------------------------------
fei::SharedPtr<LibraryWrapper>
fei::create_LibraryWrapper(MPI_Comm comm,
			       const char* libraryName)
{
  std::string libname(libraryName);

  fei::SharedPtr<LinearSystemCore> lsc;
  fei::SharedPtr<FiniteElementData> fedata;
  fei::SharedPtr<LibraryWrapper> wrapper;

  if (libname == "Aztec") {
#ifdef HAVE_AZTECOO
    lsc.reset(new Aztec_LinSysCore(comm));
#else
    std::string msg("Aztec not available.");
    throw std::runtime_error(msg);
#endif
  }

  if (libname == "FETI") {
#ifdef FEI_HAVE_FETI
    fedata.reset(new FETI_DP_FiniteElementData(comm));
#endif
  }

  if (libname == "PETSc") {
#ifdef FEI_HAVE_PETSC
    lsc.reset(new PETSc_LinSysCore(comm));
#else
    std::string msg("PETSc not available.");
    throw std::runtime_error(msg);
#endif
  }

  if (libname == "Prometheus") {
#ifdef FEI_HAVE_PROMETHEUS
    lsc.reset(new Prometheus_LinSysCore(comm));
#else
    std::string msg("Prometheus not available.");
    throw std::runtime_error(msg);
#endif
  }

  if (libname == "HYPRE") {
#ifdef FEI_HAVE_HYPRE
    lsc.reset(new HYPRE_LinSysCore(comm));
#else
    std::string msg("This factory doesn't provide Hypre instantiations.");
    throw std::runtime_error(msg);
#endif
  }

  if (libname == "TEST_LSC") {
    lsc.reset(new TEST_LinSysCore(comm));
  }

  if (lsc.get() == NULL && fedata.get() == NULL) {
    //libraryName not found
    std::string msg("create_LibraryWrapper: ");
    msg += libraryName;
    msg += " not a valid name.";
    throw std::runtime_error(msg);
  }

  if (lsc.get() != NULL) {
    wrapper.reset(new LibraryWrapper(lsc));
    return(wrapper);
  }

  if (fedata.get() != NULL) {
    wrapper.reset(new LibraryWrapper(fedata));
    return(wrapper);
  }

  return(wrapper);
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::Factory>
fei::create_fei_Factory(MPI_Comm comm,
			    const char* libraryName)
{
  std::string libname(libraryName);

  if (libname.find("Trilinos") != std::string::npos) {
    fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(comm));

    if (libname.find("Amesos") != std::string::npos) {
      fei::ParameterSet paramset;
      paramset.add(fei::Param("Trilinos_Solver", "Amesos"));
      factory->parameters(paramset);
    }
    else if (libname.find("Aztec") != std::string::npos) {

      //if libname contains "AztecOO" we'll return the Trilinos factory
      //but if libname only contains "Aztec" then we want to skip on down
      //and return an snl_fei::Factory with a Aztec LibraryWrapper...

      if (libname.find("AztecOO") != std::string::npos) {
        return(factory);
      }
    }
    else {
      //This else handles the case where libname contains "Trilinos", but
      //doesn't contain "Aztec" or "Amesos"...
      return(factory);
    }
  }

  fei::SharedPtr<LibraryWrapper> wrapper;
  try {
    wrapper = fei::create_LibraryWrapper(comm, libraryName);
  }
  catch (std::runtime_error& exc) {
    std::string msg("create_fei_Factory: ");
    msg += exc.what();
    throw std::runtime_error(msg);
  }

  if (wrapper.get() != NULL) {
    fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(comm, wrapper));
    return(factory);
  }

  fei::SharedPtr<fei::Factory> empty;
  return(empty);
}

