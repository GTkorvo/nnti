/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */

/*! \file TargetCalculator.hpp

\brief The Mesquite::TargetCalculator class is the base class. Concrete classes are 
 instantiated by the user, and often implemented by the user to give 
 mesquite a measure of the perfect mesh. 

  \author Thomas Leurent
  \date   2004-09-31
 */


#ifndef TargetCalculator_hpp
#define TargetCalculator_hpp

#include "Mesquite.hpp"
#include "MesquiteError.hpp"
#include "MsqTimer.hpp"
#include "MsqMessage.hpp"
#include "Matrix3D.hpp"
#include "PatchData.hpp"

namespace Mesquite
{
  class PatchDataParameters;

  
  /*! \class TargetCalculator
    \brief Base class that provides the interface for computing the target corner matrices
    used in the context based smoothing.

    To implement a concrete TargetCalculator, one inherits from this class and then 
    overrides the compute_target_matrices function itself to provide corner matrices for all
    elements in the patch.

    Note that an implementation is provided in TargetCalculator::compute_default_target_matrices
    for default target matrices often used in the computation of more complex,
    reference-based target matrices.

    The target calculator is set on the QualityImprover. At runtime, it associates with
    each mesh corner (i.e. a corner of an element. Most of the time, one vertex
    corresponds to many corners.) with a TargetMatrix through an MsqTag mechanism 
    available in the MsqMeshEntity class.
   */
  class TargetCalculator 
  {
  public:

    TargetCalculator() : refMesh(0), originator(0) { }
    
      //! virtual destructor ensures use of polymorphism during destruction
    virtual ~TargetCalculator()
    {};

    void initialize_default_target_matrices(Matrix3D &tri_M3D, Matrix3D &quad_M3D,
                                    Matrix3D &tet_M3D, Matrix3D &hex_M3D);

      //! \enum chooses whether the calculation is per element or an average
      //! for some cases of the \f$ \lambda_k \f$ coefficient.
    enum Lambda_type {
      REGULAR, //!< Each element has a lambda coefficient 
      AVERAGE  //!< The Lambda coefficient is the average on the mesh.
    };
      //! Computes the \f$ \lambda \f$ coefficient when it is Mesh-Based,
    double compute_L(enum Lambda_type l_type, MsqError &err);
      //! Computes the \f$ \lambda \f$ coefficient when it is element-based.
    void compute_Lk(enum Lambda_type l_type, PatchData &ref_pd, size_t elem_ind, double L_k[], int num, MsqError &err);
    
      //! \enum chooses the type of guide matrix used in the target calculator
    enum guide_type {
      Ad, //!<
      AK, //!<
      A0, //!<
      Ar, //!<
      As, //!<
      Ab, //!<
      Ac, //!<
      Ap, //!<
      Ae, //!<
      Af, //!<
      Ax  //!<
    };
      //! Computes the guide corner matrices A for a given element index in the reference patch.
    void compute_guide_matrices(enum guide_type type, PatchData &ref_pd, size_t elem_ind,
                                           Matrix3D A[], int num, MsqError &err);

    Matrix3D compute_V_3D(const Matrix3D &A, MsqError &err);
    Matrix3D compute_Q_3D(const Matrix3D &A, MsqError &err);
    Matrix3D compute_Delta_3D(const Matrix3D &A, MsqError &err);
    
      //! Compute the default "isotropic" target matrices that are often used in the computation
      //! of reference-based target matrices.
      //! The resulting corner matrices are stored in tags on the elements of the PatchData.
    void compute_default_target_matrices(PatchData &pd, MsqError &err);

      //! Compute the corner matrices for the reference mesh refMesh.
      //! The refMesh data member is set by the constructors of a concrete TargetCalculator
      //! that requires a reference mesh.
    void compute_reference_corner_matrices(PatchData &pd, MsqError &err);

    //! This function wraps compute_target_matrices and checks that the determinant of each target
    //! is positive.
    void compute_target_matrices_and_check_det(PatchData& pd, MsqError& err);
    
    /*! \brief This function provides the corner matrices for all elements on the Patch.

         Useful functionality includes: MsqMeshEntity::set_tag, MsqTag::target_matrix,
         MsqTag::scalar .
    */
    virtual void compute_target_matrices(PatchData& pd, MsqError& err) =0;

    void set_originator(PatchDataParameters* pdm, MsqError &err)
    { if (originator != 0)
        err.set_msg("Each TargetCalculator can be set on one object only.");
      else originator = pdm; }

  protected:
    MeshSet* refMesh;
    PatchDataParameters* originator; //! This is the object the TargetCalculator is attached to.

  };

  
#undef __FUNC__
#define __FUNC__ "TargetCalculator::initialize_default_target_matrices" 
  inline void TargetCalculator::initialize_default_target_matrices(Matrix3D &tri_M3D,
                                                      Matrix3D &quad_M3D,
                                                      Matrix3D &tet_M3D,
                                                      Matrix3D &hex_M3D)
  {
    const double v_tri[] = {1, 0.5, 0, 0, MSQ_SQRT_THREE/2, 0, 0, 0, 0};
    Matrix3D m1(v_tri);
    tri_M3D = m1;

    const double v_quad[] = {1, 0, 0, 0, 1, 0, 0, 0, 0};
    Matrix3D m2(v_quad);
    quad_M3D = m2;
    
    const double v_tet[] = {1, 0.5, 0.5, 0, MSQ_SQRT_THREE/2, MSQ_SQRT_THREE/6, 0, 0, MSQ_SQRT_TWO/MSQ_SQRT_THREE};
    Matrix3D m3(v_tet);
    tet_M3D = m3;

    const double v_hex[] = {1, 0, 0,  0, 1, 0,  0, 0, 1};
    Matrix3D m4(v_hex);
    hex_M3D = m4;
  }

  //!
#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_V_3D"
  inline  Matrix3D TargetCalculator::compute_V_3D(const Matrix3D &A, MsqError &err)
  {
    Vector3D a1(A[0][0], A[1][0], A[2][0]); 
    Vector3D a2(A[0][1], A[1][1], A[2][1]); 
    Vector3D a3(A[0][2], A[1][2], A[2][2]); 

    double a1_norm = A.column_length(0);
    Vector3D a1_x_a2 = a1 * a2;
    double a1_x_a2_norm = a1_x_a2.length();
    
    Matrix3D V;
    Vector3D v1, v2, v3;
    
    // note : % is the dot product
    v1 = (1/a1_norm) * a1;
    v2 = ((-(a1%a2) * a1) + (a1_norm*a1_norm)*a2) / (a1_norm * a1_x_a2_norm);
    v3 = (1/a1_x_a2_norm) * a1_x_a2;

    V.set_column(0, v1);
    V.set_column(1, v2);
    V.set_column(2, v3);

    return V;
  }
  
  //!
#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_Q_3D"
  inline  Matrix3D TargetCalculator::compute_Q_3D(const Matrix3D &A, MsqError &err)
  {
    Vector3D a1(A[0][0], A[1][0], A[2][0]); 
    Vector3D a2(A[0][1], A[1][1], A[2][1]); 
    Vector3D a3(A[0][2], A[1][2], A[2][2]); 

    double a1_norm = A.column_length(0);
    double a2_norm = A.column_length(1);
    double a3_norm = A.column_length(2);
    
    Vector3D a1_x_a2 = a1 * a2;
    double a1_x_a2_norm = a1_x_a2.length();
    Vector3D a1_x_a3 = a1 * a3;

    double nu = pow(a1_norm*a2_norm*a3_norm, 1/3);
    double det_A = det(A);
    double fac = nu * pow(det_A, 1/3); // ?? make sure this is right
    
    Matrix3D Q;

    Q[0][0] = fac * 1;
    Q[0][1] = fac * a1%a2 / (a1_norm*a2_norm);
    Q[0][2] = fac * a1%a3 / (a1_norm*a3_norm);
    Q[1][1] = fac * a1_x_a2_norm / (a1_norm*a2_norm);
    Q[1][2] = fac * a1_x_a2 % a1_x_a3 / (a1_x_a2_norm * a1_norm * a3_norm);
    Q[2][2] = fac * det_A / (a1_x_a2_norm*a3_norm); 
    
    return Q;
  }

  //!
#undef __FUNC__
#define __FUNC__ "TargetCalculator::compute_Delta_3D"
  inline  Matrix3D TargetCalculator::compute_Delta_3D(const Matrix3D &A, MsqError &err)
  {
    double a1_norm = A.column_length(0);
    double a2_norm = A.column_length(1);
    double a3_norm = A.column_length(2);
    
    double nu = pow(a1_norm*a2_norm*a3_norm, 1/3);
    double fac = 1/nu ;
    
    Matrix3D Delta;

    Delta[0][0] = fac * a1_norm;
    Delta[1][1] = fac * a2_norm;
    Delta[2][2] = fac * a3_norm;
    
    return Delta;
  }

  
} //namespace


#endif // TargetCalculator_hpp
