// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 15-Jan-03 at 08:05:56
//  LAST-MOD: 18-Feb-03 at 17:05:40 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*!
  \file   FeasibleNewton.cpp
  \brief  

  Implements the FeasibleNewton class member functions.
  
  \author Thomas Leurent
  \date   2003-01-15
*/
// DESCRIP-END.
//

#include "FeasibleNewton.hpp"
#include "MsqFreeVertexIndexIterator.hpp"

using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "FeasibleNewton::FeasibleNewton" 
FeasibleNewton::FeasibleNewton(ObjectiveFunction* of) :
  VertexMover(),
  objFunc(of)
{
  MsqError err;
  gradientLessThan=.001;
  maxIteration=6;
  maxCGiter = 1000;
  this->set_name("FeasibleNewton");
  set_patch_type(PatchData::GLOBAL_PATCH, err);
}  
  

#undef __FUNC__
#define __FUNC__ "FeasibleNewton::initialize" 
void FeasibleNewton::initialize(PatchData &pd, MsqError &err)
{
  // Cannot do anything.  Variable sizes with maximum size dependent
  // upon the entire MeshSet.
}

#undef __FUNC__
#define __FUNC__ "FeasibleNewton::initialize_mesh_iteration" 
void FeasibleNewton::initialize_mesh_iteration(PatchData &pd, MsqError &err)
{
  // Cannot do anything.  Variable sizes with maximum size dependent
  // upon the entire MeshSet.
}

#undef __FUNC__
#define __FUNC__ "FeasibleNewton::optimize_vertex_positions" 
void FeasibleNewton::optimize_vertex_positions(PatchData &pd, 
                                               MsqError &err)
{
  PRINT_INFO("\no  Performing Feasible Newton optimization.\n");
  int num_free_vertices = pd.num_free_vertices(err);
  int nv = pd.num_vertices();
  Vector3D* grad = new Vector3D[nv];
  bool fn_bool=true;// bool used for determining validity of patch
  /* Computes the value of the stopping criterion*/
  MeshSet *mesh=get_mesh_set();
  bool inner_criterion=inner_criterion_met(*mesh,err);
  int n, i;
  
  // variables used in conjugate gradient solver
  Vector3D* d = new Vector3D[nv];
  Vector3D* r = new Vector3D[nv];
  Vector3D* z = new Vector3D[nv];
  Vector3D* p = new Vector3D[nv];
  Vector3D* w = new Vector3D[nv];
  double alpha_, alpha, beta; 

  // 1.  Allocate a hessian and calculate the sparsity pattern.
  mHessian.initialize(pd, err); MSQ_CHKERR(err);
  // 2.  Calculate the gradient and Hessian for the patch
  //     (a) if not defined at current point, stop and throw an error
  fn_bool = objFunc->compute_hessian(pd, mHessian, err); MSQ_CHKERR(err);
  if (!fn_bool) { err.set_msg("invalid patch for hessian calculation"); return; }
  fn_bool = objFunc->compute_gradient(pd, grad, err); MSQ_CHKERR(err);
  if (!fn_bool) { err.set_msg("invalid patch for gradient calculation"); return; }
  // 3.  Calculate the norm of the gradient for the patch
  double grad_norm=0;
  for (n=0; n<nv; ++n) 
    grad_norm += grad[n] % grad[n]; // dot product
  grad_norm = sqrt(grad_norm);
  MSQ_DEBUG_ACTION(3,{std::cout<< "  o  gradient norm: " << grad_norm << std::endl;});

  
  // does the Feasible Newton iteration until stopping is required.
  // Terminate when: (a) too many iterations or (b) norm of the 
  // gradient of the patch is small.

  int nb_iterations = 0;
  while ( nb_iterations<maxIteration && grad_norm>gradientLessThan && !inner_criterion) {
    
    ++nb_iterations;

    // Prints out free vertices coordinates. 
    MSQ_DEBUG_ACTION(3,{
      std::cout << "\n  o Free vertices ("<< num_free_vertices
                <<")original coordinates:\n ";
      MsqVertex* toto1 = pd.get_vertex_array(err); MSQ_CHKERR(err);
      MsqFreeVertexIndexIterator ind1(&pd, err); MSQ_CHKERR(err);
      ind1.reset();
      while (ind1.next()) {
        std::cout << "\t\t\t" << toto1[ind1.value()];
      }
    });
      
    double original_value = 0.0;
    fn_bool=objFunc->evaluate(pd, original_value, err);  MSQ_CHKERR(err);
    if(!fn_bool){
      err.set_msg("Feasible Newton passed invalid patch");
    }
    MSQ_DEBUG_ACTION(3,{std::cout << "  o  original_value: " << original_value
                                  << std::endl;});
    

    // Prints out free vertices coordinates. 
    MSQ_DEBUG_ACTION(3,{
      std::cout << "  o Free vertices new coordinates: \n";
      MsqVertex* toto1 = pd.get_vertex_array(err); MSQ_CHKERR(err);
      MsqFreeVertexIndexIterator ind(&pd, err); MSQ_CHKERR(err);
      ind.reset();
      while (ind.next()) {
        std::cout << "\t\t\t" << toto1[ind.value()];
      }
    });

    // 4. Calculate a preconditioner (not needed right now)
    // 5. Calculate direction using conjugate gradients to find a
    //    zero of the Newton system of equations (H*d = -g)
    //    (a) stop if conjugate iteration limit reached
    //    (b) stop if relative residual is small
    //    (c) stop if direction of negative curvature is obtained
    double cg_tol =10e-2;
    double norm_g = length(grad, nv);
    double norm_r = norm_g;
    double rzm1; // r^T_{k-1} z_{k-1}
    double rzm2; // r^T_{k-2} z_{k-2}
    mHessian.compute_preconditionner(err); MSQ_CHKERR(err); // get M^{-1}

    for (i=0; i<nv; ++i)  d[i] = 0. ;  
    for (i=0; i<nv; ++i)  r[i] = -grad[i] ;  // r = -g because x_0 = 0 (and b=-g) 
    norm_g *= cg_tol;

    mHessian.apply_preconditionner(z, r, err); // solve Mz = r (computes z = M^-1 r)
    for (i=0; i<nv; ++i)  p[i] = z[i] ; // p_1 = z_0  
    rzm1 = inner(z,r,nv); // inner product r_{k-1}^T z_{k-1} 
    
    int cg_iter = 0;
    while ((norm_r > norm_g) && (cg_iter < maxCGiter)) {
      ++cg_iter;
      
      axpy(w, nv, mHessian, p, nv, 0,0,err); // w = A * p_k
      
      alpha_ = inner(p,w,nv); // alpha_ = p_k^T A p_k
      if (alpha_ <= 0.0) {
        printf("Direction of Negative Curvature\n");
        break; // Newton goes on with this direction of negative curvature 
      }
      
      alpha = rzm1 / alpha_;
      
      for (i=0; i<nv; ++i)  d[i] += alpha*p[i]; // x_{k+1} = d_k + alpha_{k+1} p_{k+1} 
      for (i=0; i<nv; ++i)  r[i] -= alpha*w[i]; // r_{k+1} = r_k - alpha_{k+1} A p_{k+1} 
      norm_r = length(r, nv);
      
      mHessian.apply_preconditionner(z, r, err); // solve Mz = r (computes z = M^-1 r)
      
      rzm2 = rzm1;
      rzm1 = inner(z,r,nv); // inner product r_{k-1}^T z_{k-1} 
      beta = rzm1 / rzm2;
      for (i=0; i<nv; ++i)  p[i] = z[i] + beta*p[i]; // p_k = z_{k-1} + Beta_k * p_{k-1}
    }


    // 6. Check for descent direction (inner produce of gradient and
    //    direction is negative.
    // 7. Search along the direction
    //    (a) trial = x + beta*d
    //    (b) gradient evaluation  
    //    (c) check for sufficient decrease and stop
    //    (d) otherwise, shrink beta
    // 8. Set x to trial point and calculate Hessian if needed

    bool inner_criterion=inner_criterion_met(*mesh,err);
  }

  delete[] grad;
  
}


#undef __FUNC__
#define __FUNC__ "FeasibleNewton::terminate_mesh_iteration" 
void FeasibleNewton::terminate_mesh_iteration(PatchData &pd, MsqError &err)
{
  //  std::cout << "- Executing FeasibleNewton::iteration_complete()\n";
}
  
#undef __FUNC__
#define __FUNC__ "FeasibleNewton::cleanup" 
void FeasibleNewton::cleanup()
{
  //  std::cout << "- Executing FeasibleNewton::iteration_end()\n";
}
  

