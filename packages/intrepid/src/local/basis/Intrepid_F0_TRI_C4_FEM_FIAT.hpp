#ifndef INTREPID_F0_TRI_C4_FEM_FIAT_HPP
#define INTREPID_F0_TRI_C4_FEM_FIAT_HPP
#include "Intrepid_Basis.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Intrepid_Utils.hpp"
#include "Intrepid_Tabulate.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS_types.hpp"

namespace Intrepid {

/** \class Intrepid::Basis_F0_TRI_C4_FEM_FIAT
  \brief Implementation of FIAT FEM basis functions of degree 4 for 0-forms on TRI cells. 
  Reconstruction space type is COMPLETE. Definition of the DoF set 
  for this basis, its enumeration and the associated local DoF tags are as follows,

\verbatim
  =================================================================================================
  |        |                degree-of-freedom-tag              |                                  |
  | DoF Id |---------------------------------------------------|          DoF definition          |
  |        |  subc dim  |  subc id   | subc DoFId |subc num DoF|                                  |
  |========|============|============|============|============|==================================|
  |    0   |     0      |     0      |     0      |     1      |               L_0(u)=u(0.0,0.0)  |
  |========|============|============|============|============|==================================|
  |    1   |     0      |     1      |     0      |     1      |               L_1(u)=u(1.0,0.0)  |
  |========|============|============|============|============|==================================|
  |    2   |     0      |     2      |     0      |     1      |               L_2(u)=u(0.0,1.0)  |
  |========|============|============|============|============|==================================|
  |    3   |     1      |     0      |     0      |     3      |              L_3(u)=u(0.25,0.0)  |
  |========|============|============|============|============|==================================|
  |    4   |     1      |     0      |     1      |     3      |               L_4(u)=u(0.5,0.0)  |
  |========|============|============|============|============|==================================|
  |    5   |     1      |     0      |     2      |     3      |              L_5(u)=u(0.75,0.0)  |
  |========|============|============|============|============|==================================|
  |    6   |     1      |     1      |     0      |     3      |             L_6(u)=u(0.75,0.25)  |
  |========|============|============|============|============|==================================|
  |    7   |     1      |     1      |     1      |     3      |               L_7(u)=u(0.5,0.5)  |
  |========|============|============|============|============|==================================|
  |    8   |     1      |     1      |     2      |     3      |             L_8(u)=u(0.25,0.75)  |
  |========|============|============|============|============|==================================|
  |    9   |     1      |     2      |     0      |     3      |              L_9(u)=u(0.0,0.75)  |
  |========|============|============|============|============|==================================|
  |   10   |     1      |     2      |     1      |     3      |              L_10(u)=u(0.0,0.5)  |
  |========|============|============|============|============|==================================|
  |   11   |     1      |     2      |     2      |     3      |             L_11(u)=u(0.0,0.25)  |
  |========|============|============|============|============|==================================|
  |   12   |     2      |     0      |     0      |     3      |            L_12(u)=u(0.25,0.25)  |
  |========|============|============|============|============|==================================|
  |   13   |     2      |     0      |     1      |     3      |             L_13(u)=u(0.5,0.25)  |
  |========|============|============|============|============|==================================|
  |   14   |     2      |     0      |     2      |     3      |             L_14(u)=u(0.25,0.5)  |
  |========|============|============|============|============|==================================|


\endverbatim

  The DefaultBasisFactory will select this basis
  if the following parameters are specified:
  \verbatim
  |=======================|===================================|
  |  EField               |  FIELD_FORM_0                     |
  |-----------------------|-----------------------------------|
  |  ECell                |  CELL_TRI                   |
  |-----------------------|-----------------------------------|
  |  EReconstructionSpace |  RECONSTRUCTION_SPACE_COMPLETE    |
  |-----------------------|-----------------------------------|
  |  degree               |  4                          |
  |-----------------------|-----------------------------------|
  |  EBasis               |  BASIS_FEM_FIAT                   |
  |-----------------------|-----------------------------------|
  |  ECoordinates         |  COORDINATES_CARTESIAN            |
  |=======================|===================================|
  \endverbatim

*/

template<class Scalar>
class Basis_F0_TRI_C4_FEM_FIAT: public Basis<Scalar> {
  private:
  
  /** \brief Dimension of the space spanned by the basis = number of degrees of freedom.
  */
  int numDof_;
  
  /**\brief Lookup table for the DoF's local enumeration (DoF Id) by its local DoF tag
  */
  Teuchos::Array<Teuchos::Array<Teuchos::Array<int> > > tagToEnum_;
  
  /**\brief Lookup table for the DoF's local DoF tag by its local enumeration (DoF Id)
  */
  Teuchos::Array<LocalDofTag> enumToTag_;
  
  /**\brief "true" if both lookup arrays have been set by initialize()
  */
  bool isSet_;

  /**\brief coefficients of nodal basis in terms of orthogonal polynomials */
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > vdm_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > dmats0_;
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> > dmats1_;

/* \brief Static data array that provides the Vandermonde and derivative matrices.  They need to
     be static data members because of the templating w.r.t type */
  #ifndef DOXYGEN_SHOULD_SKIP_THIS
  static Scalar *get_vdm_data() {
    static Scalar vdm_data[] = { -1.387778780781446e-17 , -2.775557561562891e-17 , 1.581627245398537e-17 , 8.888888888888893e-02 , -2.222222222222223e-02 , 8.888888888888892e-02 , 8.888888888888893e-02 , -2.222222222222225e-02 , 8.888888888888892e-02 , 8.888888888888892e-02 , -2.222222222222231e-02 , 8.888888888888904e-02 , 1.777777777777777e-01 , 1.777777777777778e-01 , 1.777777777777778e-01 , -7.142857142857151e-02 , 7.142857142857141e-02 , -5.357266658508691e-17 , -1.523809523809523e-01 , 5.979931423642378e-17 , 1.523809523809523e-01 , 2.285714285714286e-01 , 1.904761904761898e-02 , 7.619047619047616e-02 , -7.619047619047624e-02 , -1.904761904761887e-02 , -2.285714285714286e-01 , -6.095238095238094e-01 , 6.095238095238096e-01 , -2.142906663403477e-16 , -2.380952380952381e-02 , -2.380952380952382e-02 , 4.761904761904762e-02 , -1.015873015873016e-01 , -1.269841269841269e-02 , -1.015873015873016e-01 , -2.539682539682538e-02 , 6.349206349206327e-03 , 1.269841269841270e-01 , 1.269841269841269e-01 , 6.349206349206327e-03 , -2.539682539682538e-02 , -2.031746031746032e-01 , -2.031746031746032e-01 , 4.063492063492065e-01 , 7.936507936507929e-02 , 7.936507936507936e-02 , 0.000000000000000e+00 , 2.539682539682540e-01 , -6.666666666666666e-01 , 2.539682539682540e-01 , 3.809523809523810e-01 , -1.110223024625157e-16 , 0.000000000000000e+00 , -1.691768418476429e-16 , -1.110223024625157e-16 , 3.809523809523810e-01 , -3.809523809523810e-01 , -3.809523809523810e-01 , 2.220446049250313e-16 , 4.761904761904765e-02 , -4.761904761904763e-02 , 0.000000000000000e+00 , 2.285714285714286e-01 , 6.661338147750939e-17 , -2.285714285714286e-01 , -1.523809523809523e-01 , 4.000000000000000e-01 , 7.619047619047614e-02 , -7.619047619047623e-02 , -4.000000000000000e-01 , 1.523809523809524e-01 , -2.285714285714286e-01 , 2.285714285714286e-01 , 0.000000000000000e+00 , 1.587301587301587e-02 , 1.587301587301587e-02 , 4.761904761904762e-02 , 8.888888888888891e-02 , 6.666666666666667e-02 , 8.888888888888891e-02 , 3.809523809523811e-02 , -2.000000000000000e-01 , 1.904761904761905e-01 , 1.904761904761904e-01 , -2.000000000000000e-01 , 3.809523809523811e-02 , -7.619047619047623e-02 , -7.619047619047623e-02 , -2.285714285714285e-01 , -1.481481481481481e-01 , 1.481481481481482e-01 , -2.631639762074445e-17 , 2.962962962962962e-01 , -1.184237892933500e-16 , -2.962962962962962e-01 , 2.370370370370370e-01 , 0.000000000000000e+00 , 4.386066270124075e-17 , -4.386066270124075e-17 , 8.772132540248150e-17 , -2.370370370370371e-01 , 7.111111111111111e-01 , -7.111111111111110e-01 , -1.754426508049630e-16 , -1.058201058201058e-01 , -1.058201058201058e-01 , 0.000000000000000e+00 , -4.232804232804234e-02 , 2.962962962962963e-01 , -4.232804232804234e-02 , 8.465608465608462e-02 , 2.962962962962963e-01 , 0.000000000000000e+00 , -3.759485374392063e-17 , 2.962962962962963e-01 , 8.465608465608462e-02 , -8.465608465608462e-02 , -8.465608465608462e-02 , -5.925925925925926e-01 , -6.349206349206350e-02 , 6.349206349206349e-02 , 2.819614030794048e-17 , -1.269841269841270e-01 , 0.000000000000000e+00 , 1.269841269841270e-01 , -1.523809523809524e-01 , 5.639228061588096e-17 , 2.539682539682540e-01 , -2.539682539682540e-01 , -8.458842092382145e-17 , 1.523809523809524e-01 , 3.047619047619048e-01 , -3.047619047619048e-01 , 1.127845612317619e-16 , -2.116402116402115e-02 , -2.116402116402116e-02 , 8.465608465608467e-02 , -5.925925925925927e-02 , -5.925925925925927e-02 , -5.925925925925927e-02 , 6.772486772486770e-02 , -5.925925925925923e-02 , -1.693121693121694e-02 , -1.693121693121690e-02 , -5.925925925925923e-02 , 6.772486772486770e-02 , 2.201058201058201e-01 , 2.201058201058201e-01 , -2.878306878306879e-01 , 1.523809523809524e-01 , 1.523809523809524e-01 , 0.000000000000000e+00 , -6.095238095238096e-01 , 9.142857142857143e-01 , -6.095238095238096e-01 , -8.020235465369737e-17 , 8.020235465369737e-17 , 0.000000000000000e+00 , 3.561700030627496e-32 , 8.020235465369737e-17 , -8.020235465369737e-17 , 8.020235465369737e-17 , 8.020235465369737e-17 , -1.604047093073947e-16 , 1.185185185185185e-01 , -1.185185185185185e-01 , -2.631639762074445e-17 , -2.370370370370371e-01 , -1.315819881037222e-16 , 2.370370370370371e-01 , 2.370370370370370e-01 , 0.000000000000000e+00 , 4.386066270124075e-17 , -4.386066270124075e-17 , 8.772132540248150e-17 , -2.370370370370371e-01 , 7.111111111111112e-01 , -7.111111111111110e-01 , -1.754426508049630e-16 , 8.465608465608468e-02 , 8.465608465608465e-02 , 0.000000000000000e+00 , -4.232804232804233e-02 , -8.465608465608465e-02 , -4.232804232804233e-02 , -2.962962962962963e-01 , 2.962962962962963e-01 , 0.000000000000000e+00 , 1.315819881037222e-16 , 2.962962962962963e-01 , -2.962962962962963e-01 , 2.962962962962963e-01 , 2.962962962962963e-01 , -5.925925925925926e-01 , 5.079365079365080e-02 , -5.079365079365080e-02 , -2.537652627714643e-17 , 2.539682539682540e-02 , 6.344131569286609e-18 , -2.539682539682540e-02 , 2.285714285714286e-01 , -3.809523809523810e-01 , 2.539682539682540e-01 , -2.539682539682540e-01 , 3.809523809523810e-01 , -2.285714285714286e-01 , -7.619047619047613e-02 , 7.619047619047623e-02 , -1.015061051085857e-16 , 1.693121693121693e-02 , 1.693121693121693e-02 , 8.465608465608465e-02 , 1.693121693121693e-02 , 1.693121693121693e-02 , 1.693121693121693e-02 , -8.465608465608465e-02 , 1.693121693121693e-01 , -1.693121693121693e-01 , -1.693121693121693e-01 , 1.693121693121693e-01 , -8.465608465608465e-02 , -8.465608465608465e-02 , -8.465608465608465e-02 , 1.693121693121693e-01 };
    return vdm_data;
  }

  static Scalar *get_dmats0_data() {
    static Scalar dmats0_data[] = { 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 2.000000000000004e+00 , 3.655717703862945e-15 , 1.110223024625157e-15 , -2.615963001773027e-15 , 1.105782132526656e-14 , 8.881784197001252e-16 , 2.363803130039448e-16 , 3.957739485932307e-16 , 1.171902081548777e-15 , 1.776356839400250e-15 , -8.881784197001250e-16 , 3.835605607911849e-15 , 8.814965218667098e-17 , 2.221371235104169e-15 , 8.881784197001252e-16 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , -2.297368644527914e-15 , 6.000000000000004e+00 , 2.786659791809142e-15 , -5.590069302514108e-15 , 2.664535259100375e-15 , -1.598721155460227e-15 , -1.810858561866906e-15 , 3.467930675137577e-15 , -4.440892098500626e-15 , -7.919590908992791e-16 , -8.908441114415425e-15 , 3.588179136531591e-15 , 4.950129812226260e-15 , -3.996802888650564e-15 , -9.251858538543826e-18 , 1.333333333333335e+00 , 2.564615186884110e-15 , 3.333333333333335e+00 , -1.283695372222879e-16 , 1.092459456231154e-14 , 5.329070518200751e-15 , -9.664360594998842e-16 , 2.123815526736650e-15 , 1.998401444325283e-15 , 3.552713678800501e-15 , 7.771561172376129e-16 , 3.049459300862046e-15 , -2.091562519887006e-15 , 3.935740622296183e-15 , 4.440892098500626e-15 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 1.000000000000003e+00 , -2.845805601751971e-15 , -8.000000000000023e-01 , 1.000000000000001e+01 , -5.928590951498336e-15 , 2.000000000000014e-01 , -3.126398667636448e-15 , -2.164824756845977e-15 , 1.572815951552305e-15 , -1.387778780781446e-15 , 4.632546918087245e-15 , -6.430265971767091e-15 , -1.676128371899369e-15 , -4.157785227221212e-15 , 9.714451465470120e-16 , -1.195868800810527e-15 , 2.399999999999999e+00 , 2.042810365310286e-15 , 1.051820667600539e-16 , 8.400000000000009e+00 , -1.953992523340282e-15 , -1.091205198590657e-15 , 8.558740136363865e-16 , -5.773159728050813e-15 , -1.391479524196865e-15 , -3.332056852656245e-15 , 2.330206734639411e-15 , 5.080426819977399e-15 , 7.105427357600999e-15 , 8.511709855459503e-16 , 9.999999999999993e-01 , -2.681805395038991e-15 , 3.200000000000001e+00 , 1.613292832658427e-15 , 5.062616992290711e-15 , 4.200000000000006e+00 , 7.937907719836413e-17 , -7.668762744169996e-16 , -5.674473236973006e-16 , 1.776356839400250e-15 , 2.220446049250342e-16 , -3.888490726567886e-16 , 2.648216008622659e-15 , 4.316917194084153e-15 , 4.440892098500626e-15 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , -1.199315994749959e-15 , 4.000000000000008e+00 , 1.973194850143741e-15 , -1.361605960047035e-14 , -1.714285714285720e+00 , 1.855658484016334e-15 , 1.400000000000000e+01 , 6.000818036257640e-15 , 2.857142857142886e-01 , 1.718202300015125e-16 , -1.170522012650110e-14 , 5.502799861868944e-15 , 2.331211355642314e-15 , -6.328271240363392e-15 , 2.201942332173227e-15 , 8.000000000000002e-01 , 2.261380802947377e-15 , 7.999999999999980e-01 , 2.857142857142856e+00 , -1.166051382434879e-14 , -1.485714285714291e+00 , 1.262595994833550e-15 , 1.285714285714286e+01 , 5.286776307738841e-15 , 5.142857142857136e-01 , 1.640791412203535e-15 , -5.792280235697147e-15 , 2.085523112229894e-15 , -9.470202400052587e-15 , -3.552713678800501e-15 , 2.877031298653246e-16 , 1.200000000000000e+00 , -3.253633190534154e-16 , 3.668184897932611e-16 , 5.485714285714293e+00 , 2.366361075343894e-15 , -2.176269926653864e-15 , -1.102672060487993e-15 , 1.028571428571428e+01 , 3.605581441877857e-16 , -2.680321242731770e-15 , 3.781066369092739e-15 , 3.371145954981582e-15 , 2.664535259100371e-15 , -4.625929269272011e-17 , 7.999999999999999e-01 , -2.501878774698944e-15 , 2.799999999999999e+00 , -1.020513039153213e-15 , -1.065814103640156e-15 , 4.800000000000000e+00 , 7.600168156621667e-16 , -1.043844610983546e-15 , -1.504968988936321e-15 , 4.799999999999999e+00 , -1.554312234475220e-15 , -8.004352885704412e-16 , 3.148073365774497e-15 , 3.454643978291942e-15 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 };
    return dmats0_data;
  }

  static Scalar *get_dmats1_data() {
    static Scalar dmats1_data[] = { 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 0.000000000000000e+00 , 1.000000000000002e+00 , 1.827858851931473e-15 , 5.551115123125783e-16 , -1.307981500886513e-15 , 5.528910662633279e-15 , 4.440892098500626e-16 , 1.181901565019724e-16 , 1.978869742966153e-16 , 5.859510407743884e-16 , 8.881784197001252e-16 , -4.440892098500625e-16 , 1.917802803955924e-15 , 4.407482609333549e-17 , 1.110685617552084e-15 , 4.440892098500626e-16 , 3.000000000000006e+00 , 4.789687165403695e-15 , 1.776356839400250e-15 , -8.153200337091025e-16 , 2.013944566670034e-14 , 8.881784197001252e-16 , -9.776971600442707e-16 , -4.055397992727949e-16 , 1.757853122323165e-15 , 3.108624468950438e-15 , -5.329070518200750e-15 , 7.085676041417961e-15 , 1.908581317680256e-15 , 1.555700013256003e-15 , 1.776356839400250e-15 , 6.666666666666665e-01 , 3.000000000000002e+00 , -3.333333333333315e-01 , -6.803980604239456e-15 , 4.307665335545607e-15 , -1.443289932012704e-15 , -2.090910100317081e-15 , 2.257209337137497e-15 , -1.418618309243255e-15 , -3.885780586188048e-16 , -8.006934236008214e-15 , 3.556292933164406e-15 , 6.766706536662107e-16 , -5.708396718281013e-16 , -4.718447854656915e-16 , 6.666666666666657e-01 , 5.000000000000004e+00 , 1.666666666666671e+00 , -6.572947561748162e-15 , 6.350475700855894e-15 , 8.881784197001252e-16 , -2.140296567922384e-15 , 5.210102753891483e-15 , -2.109423746787797e-15 , 1.332267629550188e-15 , -9.348087504696455e-15 , 2.442462618240390e-15 , 2.191148497211589e-15 , -1.584843367652410e-15 , 1.776356839400250e-15 , -1.333333333333339e+00 , -2.959361151195280e-15 , 6.666666666666670e+00 , 9.083012120214557e-15 , 1.509903313490210e-15 , 5.329070518200751e-15 , -6.402435633658222e-15 , 3.456083156286838e-15 , 7.648203058528872e-16 , 2.664535259100376e-15 , 4.440892098500635e-15 , -1.572292614099606e-15 , -5.691691973697539e-15 , -1.239749044164735e-16 , 4.440892098500626e-15 , 5.000000000000002e-01 , 2.400000000000001e+00 , -4.000000000000001e-01 , 5.000000000000002e+00 , -6.000000000000026e-01 , 1.000000000000001e-01 , -2.015503247837709e-15 , 1.258148127127420e-15 , -2.405483220021172e-16 , -8.881784197001252e-16 , -2.532452128659149e-15 , -2.097550528330001e-15 , 2.623287389782701e-16 , -4.493627692170321e-15 , 7.910339050454240e-16 , 4.999999999999999e-01 , 1.199999999999999e+00 , 8.000000000000004e-01 , 7.000000000000005e+00 , 4.200000000000004e+00 , -7.000000000000013e-01 , -3.546495775690791e-15 , 2.262888950295884e-15 , -1.591319668629391e-15 , 1.110223024625157e-16 , 1.754390100272977e-15 , -2.607608293153362e-15 , -9.747218464460848e-16 , 9.520162436160713e-16 , 4.440892098500626e-16 , 5.000000000000009e-01 , -3.600000000000004e+00 , 1.599999999999999e+00 , 6.945986995453438e-15 , 8.400000000000009e+00 , 2.100000000000004e+00 , 1.758543624038072e-15 , -3.217539403634722e-15 , -3.614392735724119e-15 , 8.881784197001252e-16 , 8.240919520859557e-15 , -5.642185185204478e-16 , 1.676449616987502e-15 , 9.707975164493139e-15 , 1.776356839400250e-15 , 2.500000000000000e+00 , -6.504056552595708e-15 , -2.000000000000004e+00 , -8.135853102331224e-16 , -6.661338147751002e-16 , 1.050000000000000e+01 , 4.347672614741063e-15 , -6.275844041978315e-15 , -4.403884664346452e-15 , -3.552713678800501e-15 , -6.217248937900881e-15 , -2.227244763479376e-15 , 1.240539307081591e-14 , 6.536438057480606e-15 , 0.000000000000000e+00 , 4.000000000000007e-01 , 2.000000000000001e+00 , -4.000000000000004e-01 , 4.285714285714284e+00 , -8.571428571428605e-01 , 1.714285714285737e-01 , 6.999999999999997e+00 , -7.142857142857127e-01 , 1.428571428571441e-01 , -2.857142857142896e-02 , -2.558755676476370e-15 , 5.051888574510642e-16 , 1.057281834099050e-15 , -4.202194148206218e-15 , 1.710437347313132e-15 , 4.000000000000002e-01 , 1.200000000000003e+00 , 3.999999999999999e-01 , 1.428571428571425e+00 , 2.914285714285708e+00 , -7.428571428571438e-01 , 9.000000000000000e+00 , 6.428571428571431e+00 , -1.285714285714282e+00 , 2.571428571428565e-01 , -2.337841854261606e-15 , -2.214446359167745e-15 , 3.037693553488263e-16 , -2.417510636121280e-15 , -1.110223024625157e-16 , 3.999999999999987e-01 , 6.000000000000013e-01 , 1.000000000000001e+00 , -5.714285714285717e+00 , 2.742857142857146e+00 , 5.714285714285703e-01 , 4.183175938346368e-16 , 1.028571428571428e+01 , 5.142857142857144e+00 , -1.028571428571430e+00 , -2.623342001447694e-15 , 1.534939403414216e-15 , 9.024417016137114e-15 , 2.498001805406322e-17 , -2.886579864025407e-15 , 3.999999999999996e-01 , 4.200000000000004e+00 , 1.400000000000000e+00 , -6.663784658557935e-15 , -4.800000000000003e+00 , 2.400000000000003e+00 , -3.837907241356101e-15 , 1.758141325059389e-15 , 1.200000000000000e+01 , 2.399999999999999e+00 , -9.147052328536044e-15 , 6.065079228220986e-15 , 3.461223077697124e-15 , -7.931618325092893e-15 , 8.881784197001252e-16 , -1.599999999999999e+00 , 2.715816989285446e-15 , 5.600000000000003e+00 , -2.391935855732587e-15 , -7.016609515630994e-15 , -2.400000000000013e+00 , 4.223333243170086e-16 , 4.357074665788342e-15 , -8.141635513917805e-16 , 1.440000000000000e+01 , 4.440892098500570e-16 , 5.627260722491471e-15 , -6.822474683963894e-15 , -5.312417172831383e-15 , -5.773159728050814e-15 };
    return dmats1_data;
  }
#endif


  public:

  /** \brief Constructor.
  */
  Basis_F0_TRI_C4_FEM_FIAT() : numDof_(15), isSet_(false) , vdm_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TRI_C4_FEM_FIAT::get_vdm_data(),15,15,15) ) ),dmats0_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TRI_C4_FEM_FIAT::get_dmats0_data(),15,15,15) ) ),dmats1_( rcp( new Teuchos::SerialDenseMatrix<int,Scalar>(Teuchos::View,Basis_F0_TRI_C4_FEM_FIAT::get_dmats1_data(),15,15,15) ) ) { }
  
  /** \brief Initializes arrays needed for the lookup of the local enumeration (DoF Id) of a 
    degree-of-freedom by its local DoF tag and the reverse lookup of the DoF tag by DoF Id.
   */
  void initialize();
  
  
  /** \brief Returns FieldContainer with multi-indexed quantity with the values of an operator 
  applied to a set of FEM basis functions, evaluated at a set of
  points in a <strong>reference</strong> cell. The rank of the return
  FieldContainer argument depends on the number of points, number of
  basis, functions, their rank, and the rank of the operator applied to them; see 
  FieldContainer::resize(int,int,EField,EOperator,int) for summary of
  the admissible index combinations.
  In particular, the admissible range of the <var>operatorType</var> argument depends on the rank of 
  the basis functions and the space dimension.Because FEM
  reconstruction relies on COMPLETE or INCOMPLETE polynomials, i.e.,
  smooth local spaces, the admissible range of <var>operatorType</var>
  always includes VALUE, D0, D1,...,D10. If derivative order exceeds
  polynomial degree, output container is filled with 0. 
 
      \param outputValues   [out]         - FieldContainer with the computed values (see
                                            implementation for index ordering)
      \param inputPoints     [in]         - evaluation points on the reference cell  
      \param operatorType    [in]         - the operator being applied to the basis function
  */    

  void getValues(FieldContainer<Scalar>&                  outputValues,
                 const Teuchos::Array< Point<Scalar> >& inputPoints,
                 const EOperator                        operatorType) const;
  
  
  /** \brief This method is intended for FVD reconstructions and should not be used here. Its 
    invocation will throw an exception. 
  */  
  void getValues(FieldContainer<Scalar>&                  outputValues,
                 const Teuchos::Array< Point<Scalar> >& inputPoints,
                 const Cell<Scalar>&                    cell) const;

  
  int getNumLocalDof() const;
  
  
  int getLocalDofEnumeration(const LocalDofTag dofTag);

  
  LocalDofTag getLocalDofTag(int id);

  
  const Teuchos::Array<LocalDofTag> & getAllLocalDofTags();

  
  ECell getCellType() const;

  
  EBasis getBasisType() const;

  
  ECoordinates getCoordinateSystem() const;

  
  int getDegree() const;

};

} // namespace Intrepid

#include "Intrepid_F0_TRI_C4_FEM_FIATDef.hpp"

#endif

