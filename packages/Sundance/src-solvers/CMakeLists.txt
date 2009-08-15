
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#
# NONE


#
# B) Define the header and source files (and directories)
#

#
# src
#


SET(HEADERS "")
SET(SOURCES "")


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})



#
# Vectors headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../src-utils/Utilities)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Vectors)

APPEND_SET(HEADERS
  Vectors/RTOpPack_SUNDIALS_Helpers.hpp
  Vectors/RTOpPack_SUNDIALS_Ops.hpp
  Vectors/Thyra_SUNDIALS_Ops.hpp
  Vectors/TSFBCPartitionedVSBuilder.hpp
  Vectors/TSFGhostImporter.hpp
  Vectors/TSFGhostView.hpp
  Vectors/TSFIndexableVector.hpp
  Vectors/TSFLoadableBlockVector.hpp
  Vectors/TSFLoadableVector.hpp
  Vectors/TSFLinearCombination.hpp
  Vectors/TSFLinearCombinationDecl.hpp
  Vectors/TSFLinearCombinationImpl.hpp
  Vectors/TSFAccessibleVector.hpp
  Vectors/TSFProductVectorSpace.hpp
  Vectors/TSFProductVectorSpaceDecl.hpp
  Vectors/TSFProductVectorSpaceImpl.hpp
  Vectors/TSFProductVector.hpp
  Vectors/TSFProductVectorDecl.hpp
  Vectors/TSFProductVectorImpl.hpp
  Vectors/TSFRawDataAccessibleVector.hpp
  Vectors/TSFSequentialIteratorDecl.hpp
  Vectors/TSFSequentialIteratorImpl.hpp
  Vectors/TSFVector.hpp
  Vectors/TSFVectorDecl.hpp
  Vectors/TSFVectorImpl.hpp
  Vectors/TSFVectorSpace.hpp
  Vectors/TSFVectorSpaceDecl.hpp
  Vectors/TSFVectorSpaceImpl.hpp
  Vectors/TSFVectorType.hpp
  Vectors/TSFVectorTypeExtensions.hpp
  )



#
# Concrete files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Concrete)

APPEND_SET(HEADERS
  Concrete/EpetraTSFOperator.hpp
  Concrete/TSFDenseSerialVector.hpp
  Concrete/TSFEpetraGhostImporter.hpp
  Concrete/TSFEpetraGhostView.hpp
  Concrete/TSFEpetraMatrixFactory.hpp
  Concrete/TSFEpetraMatrix.hpp
  Concrete/TSFEpetraVector.hpp
  Concrete/TSFEpetraVectorSpace.hpp
  Concrete/TSFEpetraVectorType.hpp
  Concrete/TSFIfpackOperator.hpp
  Concrete/TSFMLOperator.hpp
  Concrete/TSFLAPACKGeneralMatrix.hpp
  Concrete/TSFVectorSpace2EpetraMap.hpp
  )

APPEND_SET(SOURCES
  Concrete/EpetraTSFOperator.cpp
  Concrete/TSFDenseSerialVector.cpp
  Concrete/TSFEpetraGhostImporter.cpp
  Concrete/TSFEpetraGhostView.cpp
  Concrete/TSFEpetraMatrix.cpp
  Concrete/TSFEpetraMatrixFactory.cpp
  Concrete/TSFEpetraVector.cpp
  Concrete/TSFEpetraVectorSpace.cpp
  Concrete/TSFEpetraVectorType.cpp
  Concrete/TSFIfpackOperator.cpp
  Concrete/TSFMLOperator.cpp
  Concrete/TSFLAPACKGeneralMatrix.cpp
  Concrete/TSFVectorSpace2EpetraMap.cpp
  )

#
# Operators files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Operators)

APPEND_SET(HEADERS
  Operators/Thyra_BlockOperatorBase.hpp
  Operators/Thyra_DefaultBlockOperatorDecl.hpp
  Operators/Thyra_DefaultBlockOperatorImpl.hpp
  Operators/Thyra_SetableBlockOperatorBase.hpp
  Operators/TSFBlockOperator.hpp
  Operators/TSFBlockOperatorDecl.hpp
  Operators/TSFBlockOperatorImpl.hpp
  Operators/TSFCollectivelyConfigurableMatrixFactory.hpp
  Operators/TSFEuclideanOpWithBackwardsCompatibleApply.hpp
  Operators/TSFExplicitlyTransposeableOp.hpp
  Operators/TSFExplicitlyTransposedLTIProblemFactory.hpp
  Operators/TSFHomogeneouslyBlockedLinearOp.hpp
  Operators/TSFILUFactorizableOp.hpp
  Operators/TSFIncrementallyConfigurableMatrixFactory.hpp
  Operators/TSFInverseOperator.hpp
  Operators/TSFInverseLTIOp.hpp
  Operators/TSFLinearOperator.hpp
  Operators/TSFLinearOperatorDecl.hpp
  Operators/TSFLinearOperatorImpl.hpp
  Operators/TSFLoadableBlockOperator.hpp
  Operators/TSFLoadableMatrix.hpp
  Operators/TSFLTIProblemFactoryBase.hpp
  Operators/TSFDefaultLTIProblemFactory.hpp
  Operators/TSFMatrixFactory.hpp
  Operators/TSFMultiVectorOperator.hpp
  Operators/TSFNonlinearOperator.hpp
  Operators/TSFNonlinearOperatorBase.hpp
  Operators/TSFNonmemberOpHelpers.hpp
  Operators/TSFOpWithBackwardsCompatibleApply.hpp
  Operators/TSFPartitionedMatrixFactory.hpp
  Operators/TSFPartitionedToMonolithicConverter.hpp
  Operators/TSFRowAccessibleOp.hpp
  Operators/TSFSimplifiedLinearOpBase.hpp 
  )

APPEND_SET(SOURCES
  Operators/TSFPartitionedMatrixFactory.cpp
  )

#
# LinearSolvers files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/LinearSolvers)

APPEND_SET(HEADERS
  LinearSolvers/TSFAmesosSolver.hpp
  LinearSolvers/TSFAnasaziEigensolver.hpp
  LinearSolvers/TSFAztecSolver.hpp
  LinearSolvers/TSFAztecSolver.hpp
  LinearSolvers/TSFBelosSolver.hpp
  LinearSolvers/TSFBICGSTABSolver.hpp
  LinearSolvers/TSFBlockTriangularSolver.hpp
  LinearSolvers/TSFEigensolver.hpp
  LinearSolvers/TSFEigensolverBase.hpp
  LinearSolvers/TSFGMRESSolver.hpp
  LinearSolvers/TSFIterativeSolver.hpp
  LinearSolvers/TSFKrylovSolver.hpp
  LinearSolvers/TSFLinearSolver.hpp
  LinearSolvers/TSFLinearSolverDecl.hpp
  LinearSolvers/TSFLinearSolverImpl.hpp
  LinearSolvers/TSFLinearSolverBase.hpp
  LinearSolvers/TSFLinearSolverBaseDecl.hpp
  LinearSolvers/TSFLinearSolverBaseImpl.hpp
  LinearSolvers/TSFLinearSolverBuilder.hpp
  LinearSolvers/TSFSolverState.hpp 
  )

APPEND_SET(SOURCES
  LinearSolvers/TSFAmesosSolver.cpp
  LinearSolvers/TSFAztecSolver.cpp
  LinearSolvers/TSFBelosSolver.cpp
  LinearSolvers/TSFLinearSolverBuilder.cpp
  )

#
# Non-LinearSolvers files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/NonlinearSolvers)

APPEND_SET(HEADERS
  NonlinearSolvers/NOX_StatusTest_SafeCombo.H
  NonlinearSolvers/NOX_TSF_Group.H
  NonlinearSolvers/NOX_TSF.H
  NonlinearSolvers/NOX_TSF_StatusTestBuilder.H
  NonlinearSolvers/NOX_TSF_Vector.H
  NonlinearSolvers/TSFNOXSolver.H 
  )

APPEND_SET(SOURCES
  NonlinearSolvers/NOX_StatusTest_SafeCombo.C
  NonlinearSolvers/NOX_TSF_Group.C
  NonlinearSolvers/NOX_TSF_StatusTestBuilder.C
  NonlinearSolvers/NOX_TSF_Vector.C
  NonlinearSolvers/TSFNOXSolver.C       
  )

#
# testing files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Testing)

APPEND_SET(HEADERS
  Testing/Thyra_TestSpecifier.hpp
  Testing/TSFCompoundTester.hpp
  Testing/TSFHeatOperator1D.hpp
  Testing/TSFLinearCombinationTester.hpp
  Testing/TSFMatrixLaplacian1D.hpp
  Testing/TSFOperatorBuilder.hpp
  Testing/TSFPoissonBoltzmannJacobian.hpp
  Testing/TSFPoissonBoltzmannOp.hpp
  Testing/TSFRandomBlockMatrix.hpp
  Testing/TSFRandomSparseMatrix.hpp
  Testing/TSFTesterBase.hpp
  Testing/TSFVectorTester.hpp 
  )

APPEND_SET(SOURCES
  Testing/TSFHeatOperator1D.cpp
  Testing/TSFMatrixLaplacian1D.cpp
  Testing/TSFPoissonBoltzmannJacobian.cpp
  Testing/TSFPoissonBoltzmannOp.cpp 
  )

#
# Adaptors files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Adapters)

APPEND_SET(HEADERS
  Adapters/TSF_NVector.hpp
  )

APPEND_SET(SOURCES
  Adapters/TSF_NVector.cpp
  )

#
# Preconditioners files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/Preconditioners)

APPEND_SET(HEADERS
  Preconditioners/TSFGenericLeftPreconditioner.hpp
  Preconditioners/TSFGenericRightPreconditioner.hpp
  Preconditioners/TSFILUKPreconditionerFactory.hpp
  Preconditioners/TSFParameterListPreconditionerFactory.hpp
  Preconditioners/TSFPreconditioner.hpp
  Preconditioners/TSFPreconditionerBase.hpp
  Preconditioners/TSFPreconditionerFactory.hpp
  Preconditioners/TSFPreconditionerFactoryBase.hpp 
  )

APPEND_SET(SOURCES
  Preconditioners/TSFParameterListPreconditionerFactory.cpp
  )



IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  TSFExtended
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
