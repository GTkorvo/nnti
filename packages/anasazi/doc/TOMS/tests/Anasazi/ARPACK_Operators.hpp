#ifndef ARPACK_OPERATORS_HPP
#define ARPACK_OPERATORS_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziOperator.hpp"
#include "MyMultiVec.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_RefCountPtr.hpp"


/******************************************************************************/
/*! \class OPA< ScalarType >
  \brief Implementation of Anasazi::Operator< ScalarType > for the identity
  operator.
*/
template <class ScalarType>
class OPA : public Anasazi::Operator<ScalarType>
{
  
public:
  
  OPA()  {}
  ~OPA() {}
  
  void Apply(const Anasazi::MultiVec<ScalarType>& X, 
        Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    const ScalarType ZERO = ScalarTraits<ScalarType>::zero();
  
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    TEST_FOR_EXCEPTION(MyX == 0,Anasazi::OperatorError,"Casting failure.");
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    TEST_FOR_EXCEPTION(MyY == 0,Anasazi::OperatorError,"Casting failure.");
      
    TEST_FOR_EXCEPTION(X.GetNumberVecs() != Y.GetNumberVecs(),Anasazi::OperatorError,"Invalid input multivectors.");
    TEST_FOR_EXCEPTION(X.GetVecLength() != Y.GetVecLength(),Anasazi::OperatorError,"Invalid input multivectors.");
    
    // MyY = ONE*MyX
    MyY->MvAddMv( ONE, *MyX, ZERO, *MyX );
  }
};

/******************************************************************************/
/*! \class OPM< ScalarType >
  \brief Implementation of Anasazi::Operator< ScalarType > for the
  application of the central difference discretization of 2-D Laplacian.
*/
template <class ScalarType>
class OPM : public Anasazi::Operator<ScalarType>
{
private:
  void tv(int nx, const ScalarType *x, ScalarType *y) const {
    ScalarType dd = 4.0,
               dl = -ScalarTraits<ScalarType>::one(),
               du = -ScalarTraits<ScalarType>::one();
    
    // Compute the matrix vector multiplication y<---T*x
    // where T is a nx by nx tridiagonal matrix with DD on the 
    // diagonal, DL on the subdiagonal, and DU on the superdiagonal.
    int j;
    j = 0;
    y[j] = dd*x[j] + du*x[j+1];
    for (j=1; j<nx-1; j++) {
      y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
    }
    j = nx-1;
    y[j] = dl*x[j-1] + dd*x[j];
  }
public:
  
  OPM()  {}
  ~OPM() {}
  
  void Apply(const Anasazi::MultiVec<ScalarType>& X, 
        Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    BLAS<int,ScalarType> blas;
    int n, nx;
  
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    TEST_FOR_EXCEPTION(MyX == 0,Anasazi::OperatorError,"Casting failure.");
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    TEST_FOR_EXCEPTION(MyY == 0,Anasazi::OperatorError,"Casting failure.");
      
    TEST_FOR_EXCEPTION(X.GetNumberVecs() != Y.GetNumberVecs(),Anasazi::OperatorError,"Invalid input multivectors.");
    TEST_FOR_EXCEPTION(X.GetVecLength() != Y.GetVecLength(),Anasazi::OperatorError,"Invalid input multivectors.");
    
    int nvecs = X.GetNumberVecs();
    
    // deduce the size of the operator from the vector length...
    n = X.GetVecLength();
    // ... and the number of interior points in the discretization from that
    nx = ScalarTraits<int>::squareroot(n);
    TEST_FOR_EXCEPTION(nx*nx != n,Anasazi::OperatorError,"Invalid input.");
    
    // The rest is stolen from the ARPACK codes (see notice above)
    //
    // The matrix used is the 2 dimensional discrete Laplacian on unit
    // square with zero Dirichlet boundary condition.
    //
    // Computes y <--- OP*x, where OP is the nx*nx by nx*nx block 
    // tridiagonal matrix
    //
    //              | T -I          | 
    //              |-I  T -I       |
    //         OP = |   -I  T       |
    //              |        ...  -I|
    //              |           -I T|
    //
    // The subroutine tv is called to computed y<---T*x.
    int p, j, lo;
    for (p=0; p<nvecs; p++) {
      lo = 0;
      tv(nx,&(*MyX)[p][lo],&(*MyY)[p][lo]);
      blas.AXPY(nx,-ONE,&(*MyX)[p][lo+nx],1,&(*MyY)[p][lo],1);
    }
    
    for (j=1; j<nx-1; j++) {
      lo = j*nx;
      for (p=0; p<nvecs; p++) {
        tv(nx,&(*MyX)[p][lo],&(*MyY)[p][lo]);
        blas.AXPY(nx,-ONE,&(*MyX)[p][lo-nx],1,&(*MyY)[p][lo],1);
        blas.AXPY(nx,-ONE,&(*MyX)[p][lo+nx],1,&(*MyY)[p][lo],1);
      }
    }
    
    for (p=0; p<nvecs; p++) {
      lo = (nx-1)*nx;
      tv(nx,&(*MyX)[p][lo],&(*MyY)[p][lo]);
      blas.AXPY(nx,-ONE,&(*MyX)[p][lo-nx],1,&(*MyY)[p][lo],1);
    }
    
    // scale the vector by (1/h^2), where h is the mesh size
    ScalarType h2 = ONE/(ScalarType)((nx+1)*(nx+1));
    for (p=0; p<nvecs; p++) {
      blas.SCAL(n,ONE/h2,&(*MyY)[p][0],1);
    }
    
  }
  
};
#endif //ARPACK_OPERATORS_HPP
