/*!
  \file   ConjugateGradient.cpp
  \brief  

  The Conjugate Gradient class is a concrete vertex mover
  which performs conjugate gradient minimizaiton.

  \author Michael Brewer
  \date   2002-06-19
*/

#include "ConjugateGradient.hpp"
#include <math.h>
#include "MsqMessage.hpp"
#include "MsqTimer.hpp"
#include "MsqFreeVertexIndexIterator.hpp"

using namespace Mesquite;

#undef __FUNC__
#define __FUNC__ "ConjugateGradient::ConjugateGradient" 
ConjugateGradient::ConjugateGradient(ObjectiveFunction* objective,
                                     MsqError &err) :
  VertexMover(),
  fGrad(NULL),
  pGrad(NULL),
  pMemento(NULL),
  fNewGrad(NULL)
{
  this->set_name("ConjugateGradient");
  this->set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err, 1); MSQ_CHKERR(err);
  objFunc=objective;
    //Michael:: default to global?
  set_debugging_level(0);
    //set the default inner termination criterion
  TerminationCriterion* default_crit=get_inner_termination_criterion();
  if(default_crit==NULL){
    err.set_msg("QualityImprover did not create a default inner termination criterion.");
  }
  else{
    default_crit->add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,5,err);
  }
  
}  


#undef __FUNC__
#define __FUNC__ "ConjugateGradient::~ConjugateGradient" 
ConjugateGradient::~ConjugateGradient()
{
  // Checks that cleanup() has been called. 
  assert(fGrad==NULL);
  assert(pGrad==NULL);
  assert(pMemento==NULL);
  assert(fNewGrad==NULL);
}
  
#undef __FUNC__
#define __FUNC__ "ConjugateGradient::initialize" 
void ConjugateGradient::initialize(PatchData &pd, MsqError &err)
{
  PRINT_INFO("\no   Performing Conjugate Gradient optimization.\n");
  arraySize=5;
  fGrad = new Vector3D[ arraySize ];
  pGrad = new Vector3D[ arraySize ];
  fNewGrad = new Vector3D[ arraySize ];
  //mCoord = new Vector3D[ arraySize ];
  pMemento=pd.create_vertices_memento(err);
}


/*! \fn ConjugateGradient::set_patch_type(PatchData::PatchType type, MsqError &err)

    ConjugateGradient supports GLOBAL_PATCH and ELEMENTS_ON_VERTEX_PATCH
*/
void ConjugateGradient::set_patch_type(PatchData::PatchType type, MsqError &err,
				       int patch_param1, int patch_param2)
{
  if (type == PatchData::GLOBAL_PATCH || type == PatchData::ELEMENTS_ON_VERTEX_PATCH) {
    PatchDataUser::set_patch_type(type, err, patch_param1, patch_param2);
  } else {
    err.set_msg("Type not supported by ConjugateGradient algorythm.");
  }
}

#undef __FUNC__
#define __FUNC__ "ConjugateGradient::initialize_mesh_iteration" 
void ConjugateGradient::initialize_mesh_iteration(PatchData &/*pd*/,
                                                  MsqError &/*err*/)
{
  
}

#undef __FUNC__
#define __FUNC__ "ConjugateGradient::optimize_vertex_position"
/*!Performs Conjugate gradient minimization on the PatchData, pd.*/
void ConjugateGradient::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err){
  // pd.reorder();

  FUNCTION_TIMER_START(__FUNC__);

  Timer c_timer;
  int num_local_vertices = pd.num_vertices();
  if(num_local_vertices>arraySize){
    delete []fGrad;
    delete []pGrad;
    delete []fNewGrad;
      //Increase number to try to avoid reallocating
    arraySize=num_local_vertices + 5;
    fGrad = new Vector3D[ arraySize ];
    pGrad = new Vector3D[ arraySize ];
    fNewGrad = new Vector3D[ arraySize ];
  }
    //zero out arrays
  int zero_loop=0;
  while(zero_loop<arraySize){
    fGrad[zero_loop].set(0,0,0);
    pGrad[zero_loop].set(0,0,0);
    fNewGrad[zero_loop].set(0,0,0);
    ++zero_loop;
  }
  
  // gets the array of vertices for the patch  
  MsqVertex* vertices=pd.get_vertex_array(err);
  int ind;
    //Michael cull list:  possibly set soft_fixed flags here

   int num_vert=pd.num_vertices();
   MsqFreeVertexIndexIterator free_iter(&pd, err);
  
   MSQ_DEBUG_ACTION(1,{ 
      if(pd.num_free_vertices(err)<1) 
	PRINT_INFO("\nEmpty free vertex list in ConjugateGradient\n"); 
   });
      
   double f=0;
   //Michael, this isn't equivalent to CUBIT because we only want to check
   //the objective function value of the 'bad' elements
   //if invalid initial patch set an error.
  if( ! objFunc->compute_gradient(pd, fGrad, f, err, num_vert) ){
    MSQ_CHKERR(err);
    err.set_msg("Conjugate Gradient not able to get valid gradient and function values on intial patch.");
  }
  double grad_norm=MSQ_MAX_CAP;
  
  if(conjGradDebug>0){
    PRINT_INFO("\nCG's DEGUB LEVEL = %i \n",conjGradDebug);
    grad_norm=Linf(fGrad,num_vert);
    PRINT_INFO("\nCG's FIRST VALUE = %f,grad_norm = %f",f,grad_norm);
    PRINT_INFO("\n   TIME %f",c_timer.since_birth());
    grad_norm=MSQ_MAX_CAP;
  }
  MSQ_CHKERR(err);
  ind=0;
    //Initializing pGrad (search direction).  
  while(ind<num_vert){    
    pGrad[ind]=(-fGrad[ind]);
    ++ind;
  }
  int j=0; // total nb of step size changes ... not used much
  int i=0; // iteration counter
  int m=0; // 
  double alp=MSQ_MAX_CAP; // alp: scale factor of search direction
    //we know inner_criterion is false because it was checked in
    //loop_over_mesh before being sent here.
  bool inner_criterion=false;
  TerminationCriterion* term_crit=get_inner_termination_criterion();
  
    //while ((i<maxIteration && alp>stepBound && grad_norm>normGradientBound)
    //     && !inner_criterion){
  while(!inner_criterion){
    ++i;
      //std::cout<<"\Michael delete i = "<<i;
    int k=0;
    alp=get_step(pd,f,k,err);
    j+=k;
    if(conjGradDebug>2){
      PRINT_INFO("\n  Alp initial, alp = %20.18f",alp);
    }
    MSQ_CHKERR(err);
     // if alp == 0, revert to steepest descent search direction
    if(alp==0){
      free_iter.reset();
      while (free_iter.next()) {
        m=free_iter.value();
        pGrad[m]=(-fGrad[m]);
      }
      alp=get_step(pd,f,k,err);
      j+=k;
      if(conjGradDebug>1){
        PRINT_INFO("\n CG's search direction reset.");
        if(conjGradDebug>2)
          PRINT_INFO("\n  Alp was zero, alp = %20.18f",alp);
      }
      
    }
    if(alp!=0){
      free_iter.reset();
      while (free_iter.next()) {
        m=free_iter.value();
        vertices[m] += (alp * pGrad[m]);
          //Added move_to_ownever
        pd.snap_vertex_to_domain(m,err);
      }
      if (! objFunc->compute_gradient(pd, fNewGrad, f, err, num_vert)){
        MSQ_CHKERR(err);
        err.set_msg("Error inside Conjugate Gradient, vertices moved making function value invalid.");
      }
      
      if(conjGradDebug>0){
        grad_norm=Linf(fNewGrad,num_vert);
        PRINT_INFO("\nCG's VALUE = %f,  iter. = %i,  grad_norm = %f,  alp = %f",f,i,grad_norm,alp);
        PRINT_INFO("\n   TIME %f",c_timer.since_birth());
      }
      double s11=0;
      double s12=0;
      double s22=0;
      free_iter.reset();
      while (free_iter.next()) {
        m=free_iter.value();
        s11+=fGrad[m]%fGrad[m];
        s12+=fGrad[m]%fNewGrad[m];
        s22+=fNewGrad[m]%fNewGrad[m];
      }
      
        // Steepest Descent (takes 2-3 times as long as P-R)
        //double bet=0;
      
        // Fletcher-Reeves (takes twice as long as P-R)
        //double bet = s22/s11;

        // Polack-Ribiere        
      double bet = (s22-s12)/s11;
      free_iter.reset();
      while (free_iter.next()) {
        m=free_iter.value();
        pGrad[m]=(-fNewGrad[m]+(bet*pGrad[m]));
        fGrad[m]=fNewGrad[m];
      }
      if(conjGradDebug>2){
        PRINT_INFO(" \nSEARCH DIRECTION INFINITY NORM = %e",
                   Linf(fNewGrad,num_vert));
      }
      
    }//end if on alp == 0
      //Removing the following line of code (4/2/03) as it should not
      //be needed with the new version of Termination Criterion.
      //Update mesh before checking criterion
      //pd.update_mesh(err);
    inner_criterion=term_crit->terminate_with_function_and_gradient(pd,objFunc,
                                                                    f,fNewGrad,
                                                                    err);
    MSQ_CHKERR(err);
  }//end while
  if(conjGradDebug>0){
    PRINT_INFO("\nConjugate Gradient complete i=%i ",i);
    PRINT_INFO("\n-  FINAL value = %f, alp=%4.2e grad_norm=%4.2e",f,alp,grad_norm);
    PRINT_INFO("\n   FINAL TIME %f",c_timer.since_birth());
  }
  FUNCTION_TIMER_END();
}


#undef __FUNC__
#define __FUNC__ "ConjugateGradient::terminate_mesh_iteration" 
void ConjugateGradient::terminate_mesh_iteration(PatchData &/*pd*/,
                                                 MsqError &/*err*/)
{
    //  cout << "- Executing ConjugateGradient::iteration_complete()\n";
}


#undef __FUNC__
#define __FUNC__ "ConjugateGradient::cleanup" 
void ConjugateGradient::cleanup()
{
    //  cout << "- Executing ConjugateGradient::iteration_end()\n";
  delete []fGrad; fGrad = NULL;
  delete []pGrad; pGrad = NULL;
  delete []fNewGrad; fNewGrad = NULL;
  pMemento->~PatchDataVerticesMemento(); pMemento = NULL;
}

//!Computes a distance to move vertices given an initial position and search direction (stored in data member pGrad).
/*!Returns alp, the double which scales the search direction vector
  which when added to the old nodal positions yields the new nodal
  positions.*/
/*!\todo Michael NOTE:  ConjugateGradient::get_step's int &j is only
  to remain consisitent with CUBIT for an initial test.  It can be
  removed.*/

double ConjugateGradient::get_step(PatchData &pd,double f0,int &j,
                                   MsqError &err)
{
  size_t num_vertices=pd.num_vertices();
    //initial guess for alp
  double alp=1.0;
  int jmax=100;
  double rho=0.5;
    //feasible=false implies the mesh is not in the feasible region
  bool feasible=false;
  int found=0;
    //f and fnew hold the objective function value
  double f=0;
  double fnew=0;
    //Counter to avoid infinitly scaling alp
  j=0;
  //save memento
  pd.recreate_vertices_memento(pMemento, err);
    //if we must check feasiblility
    //while step takes mesh into infeasible region and ...
  while (j<jmax && !feasible && alp>MSQ_MIN) {
    ++j;
    pd.set_free_vertices_constrained(pMemento,pGrad,num_vertices,alp,err);
    feasible=objFunc->evaluate(pd,f,err);      
    MSQ_CHKERR(err);
      //if not feasible, try a smaller alp (take smaller step)
    if(!feasible){
      alp*=rho;
    }
  }//end while ...
  
    //if above while ended due to j>=jmax, no valid step was found.
  if(j>=jmax){
    PRINT_INFO("\nFeasible Point Not Found");
    return 0.0;
  }
    //PRINT_INFO("\nOriginal f %f, first new f = %f, alp = %f",f0,f,alp);
    //if new f is larger than original, our step was too large
  if(f>=f0){
    j=0;
    while (j<jmax && found == 0){
      ++j;
      alp *= rho;
      pd.set_free_vertices_constrained(pMemento,pGrad,num_vertices,alp,err);
        //Get new obj value
        //if patch is now invalid, then the feasible region is  convex or
        //we have an error.  For now, we assume an error.
      if(! objFunc->evaluate(pd,f,err) ){
        err.set_msg("Non-convex feasiblility region found.");
      }
      pd.set_to_vertices_memento(pMemento,err);MSQ_CHKERR(err);
        //if our step has now improved the objective function value
      if(f<f0){
        found=1;
      }
    }//   end while j less than jmax
      //PRINT_INFO("\nj = %d found = %d f = %20.18f f0 = %20.18f\n",j,found,f,f0);
      //if above ended because of j>=jmax, take no step
    if(found==0){
        //PRINT_INFO("alp = %10.8f, but returning zero\n",alp);
      alp=0.0; 
      return alp;
    }

    j=0;
      //while shrinking the step improves the objFunc value further,
      //scale alp down.  Return alp, when scaling once more would
      //no longer improve the objFunc value.  
    while(j<jmax){
      ++j;
      alp*=rho;
      //step alp in search direction from original positions
      pd.set_free_vertices_constrained(pMemento,pGrad,num_vertices,alp,err);MSQ_CHKERR(err);

        //get new objective function value
      if (! objFunc->evaluate(pd,fnew,err))
        err.set_msg("Non-convex feasiblility region found while computing new f.");
      if(fnew<f){
        f=fnew;
      }
      else{
	//Reset the vertices to original position
	pd.set_to_vertices_memento(pMemento,err);MSQ_CHKERR(err);
	alp/=rho;
	return alp;
      }
    }
    //Reset the vertices to original position and return alp
    pd.set_to_vertices_memento(pMemento,err);MSQ_CHKERR(err);
    return alp;
  }
    //else our new f was already smaller than our original
  else{
    j=0;
      //check to see how large of step we can take
    while (j<jmax && found == 0) {
      ++j;
        //scale alp up (rho must be less than 1)
      alp /= rho;
      //step alp in search direction from original positions
      pd.set_free_vertices_constrained(pMemento,pGrad,num_vertices,alp,err);MSQ_CHKERR(err);

      feasible = objFunc->evaluate(pd,fnew, err);
      if ( ! feasible ){
         alp *= rho;
	 
	 //Reset the vertices to original position and return alp
	 pd.set_to_vertices_memento(pMemento,err);MSQ_CHKERR(err);
         return alp;
      }
      if (fnew<f) { 
          f = fnew; 
       }
       else {
         found=1;
         alp *= rho;
       }
    }

    //Reset the vertices to original position and return alp
    pd.set_to_vertices_memento(pMemento,err);MSQ_CHKERR(err);
    return alp;
  }
}

/*!Quadratic one-dimensional line search.*/
/*
double ConjugateGradient::get_step(PatchData &pd,double f0,int &j,
                                   MsqError &err){
  const double CGOLD = 0.3819660;
  const double ZEPS = 1.0e-10;
  int n=pd.num_free_vertices();
  MsqVertex* vertices=pd.get_vertex_array(err);
  double a,b,d,etemp,fb,fu,fv,fw,fx,p,q,r,tol,tol1,tol2,u,v,w,x,xm;
  double e=0.0;
  d=0.0;
  tol=.001;
  int iter, maxiter;
  maxiter=100;
  a=0;
  b=.125;
  int m=0;
  fb=f0-1.0;
  iter=0;
  //find b such that a b 'should' bracket the min
  while (fb<=f0 && iter<maxiter){
    ++iter;
    b*=2.0;
    for(m=0;m<n;++m){
      mCoord[m]=mCoord[m] + (b*pGrad[m]);
      vertices[m]=(mCoord[m]);
    }
    fb=objFunc->evaluate(pd,err);
  }
  iter=0;
  x=w=v=(b/2.0);
  for(m=0;m<n;++m){
    mCoord[m]=mCoord[m] + (x*pGrad[m]);
    vertices[m]=(mCoord[m]);
  }
  fw=fv=fx=objFunc->evaluate(pd,err);
  for(iter=0;iter<maxiter;++iter){
      //PRINT_INFO("a=%f,b=%f,x=%f,iter=%i\n",a,b,x,iter);
    xm=(a+b)*.5;
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if(fabs(x-xm)<= (tol2-0.5*(b-a))){
      return x;
    }
    if(fabs(e)>tol1){
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if(q>0.0)
        p=-p;
      q=fabs(q);
      etemp=e;
      e=d;
      if(fabs(p)>=fabs(0.5*q*etemp)||(p<=q*(a-x))||(p>=q*(b-x))){
        d=CGOLD*(e=(x>=xm?a-x:b-x));
      }
      else{
        d=p/q;
        u=x+d;
        if(u-a<tol2||b-u<tol2)
        {
          if(tol1<0.0)
            d=x-xm;
          else
            d=xm-x;
        }
      }
    }
    
    else{
      d=CGOLD*(e=(x>=xm?a-x:b-x));
    }
    if(tol<0.0)
      u=(fabs(d)>=tol1?x+d:x-d);
    else
      u=(fabs(d)>=tol1?x+d:x+d);
    for(m=0;m<n;++m){
      mCoord[m]=mCoord[m] + (u*pGrad[m]);
      vertices[m]=(mCoord[m]);
    }
    fu=objFunc->evaluate(pd,err);
    if(fu<fx){
      if(u>=x)
        a=x;
      else
        b=x;
      v=w;
      w=x;
      x=u;
      fv=fw;
      fw=fx;
      fx=fu;
    }
    else{
      if(u<x)
        a=u;
      else
        b=u;
      if(fu<=fw||w==x){
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      }
      else if (fu<=fv||v==x||v==w){
        v=u;
        fv=fu;
      }
    }
  }
  for(m=0;m<n;++m){   
    vertices[m]=(mCoord[m]);
  }
    //PRINT_WARNING("TOO MANY ITERATIONS IN QUADRATIC LINE SEARCH");
  return x;
}
*/
          



