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

//       //!
//     bool get_next_reference_patch(MsqError &err)
//     { refMesh->get_next_patch(refPatch, *originator, err); }
//       //!
//     void reset_reference_meshset(MsqError &err) { refMesh->reset(err); }

      //! \enum chooses the calculation for the \f$ \lambda_k \f$ coefficient.
    enum lambda_type {
      L00, //!< Returns the scalar 1 . 
      L11, //!< 
      L12, //!< 
      L13, //!< 
      L21, //!< 
      L22, //!< 
      L31, //!< 
      L32, //!< 
      L41  //!< 
    };
      //! Computes the \f$ \lambda \f$ coefficient when it is Mesh-Based,
    double compute_L(enum lambda_type l_type, MsqError &err);
      //! Computes the \f$ \lambda \f$ coefficient when it is element-based.
    void compute_Lk(enum lambda_type l_type, PatchData &ref_pd, size_t elem_ind, double L_k[], int num, MsqError &err);
    
      //! \enum chooses the calculation for the \f$ D_k \f$ diagonal matrix.
    enum D_type {
      D00, //!< Identity matrix
      D11, //!<
      D21, //!<
      D31, //!<
      D41, //!<
      D42, //!<
      D43, //!<
      D51, //!<
      D52, //!<
      D53  //!<
    };
      //! Computes the \f$ D \f$ diagonal matrix when it is Mesh-Based,
    Matrix3D compute_D(enum D_type l_type, MsqError &err);
      //! Computes the \f$ D \f$ diagonal matrix when it is element-based ( \f$ D_k \f$ ).
    void compute_Dk(enum D_type l_type, PatchData &ref_pd, size_t elem_ind, Matrix3D D_k[], int num, MsqError &err);
    
      //! chooses the calculation for the \f$ R_k \f$ matrix.
    enum R_type {
      R00, //!<
      R11, //!<
      R21, //!<
      R31, //!<
      R32, //!<
      R33, //!<
      R41, //!<
      R42, //!<
      R43, //!<
      R44  //!<
    };
      //! Computes the \f$ R \f$ matrix when it is Mesh-Based,
    Matrix3D compute_R(enum R_type l_type, MsqError &err);
      //! Computes the \f$ R \f$ matrix when it is element-based ( \f$ R_k \f$ ).
    void compute_Rk(enum R_type l_type, PatchData &ref_pd, size_t elem_ind, Matrix3D R_k[], int num, MsqError &err);
    
      //! chooses the calculation for the \f$ R_k \f$ matrix.
    enum W_type {
      W00, //!< W matrices are set to the default corner matrices for ideal elements. 
      W11, //!<
      W21, //!< W matrices are set to the corners of the reference mesh.
      W31, //!<
      W41, //!<
      W42  //!<
    };
      //! Computes the \f$ W \f$ matrix when it is element-based ( \f$ W_k \f$ ).
    void compute_Wk(enum W_type l_type, PatchData &ref_pd, size_t elem_ind, Matrix3D W_k[], int num, MsqError &err);
    
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
    void compute_target_matrices(PatchData& pd, MsqError& err);

    void set_originator(PatchDataParameters* pdm, MsqError &err)
    { if (originator != 0)
        err.set_msg("Each TargetCalculator can be set on one object only.");
      else originator = pdm; }

  protected:
    MeshSet* refMesh;
//    PatchData refPatch;

    enum lambda_type mLambda; 
    enum D_type mD;
    enum R_type mR;
    enum W_type mW;
    
  private:
    PatchDataParameters* originator; //! This is the object the TargetCalculator is attached to.

  };

  
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


  
} //namespace


#endif // TargetCalculator_hpp
