//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/* This file provides a simple MEX interface to use ML in MATLAB.  This code is
   was originally based off of AdaptiveSA.cpp by Marzio Sala (although the two
   are radically different at this point), and is intended to be used with its
   MATLAB front-end ml.m. 

   By: Chris Siefert <csiefer@sandia.gov>
   Version History
   08/30/2006 - Added ML_epetra interface functionality.
   08/15/2006 - Added operator complexity handling.
   08/08/2006 - Moved a few variable declarations to allow sun's less forgiving
                compiler to compile the code.
   08/07/2006 - Added status checking functionality.
   08/03/2006 - Added functionality to allow ML to have multiple systems ready
                to be solved at once.  This entailed adding the system handle stuff.
   05/22/2006 - Teuchos-friendly version 
   05/16/2006 - Initial Version 
*/

#include <stdio.h>
#include <string.h>

/* Needed for ML - Horked from the examples/AdaptiveSA.cpp.  Thanks, Marzio!*/
#include "ml_config.h"
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

/* MLAPI Headers */
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_MultiLevelAdaptiveSA.h"
#include "MLAPI_DistributedMatrix.h"
#include "MLAPI_Krylov.h"


/* ML_Epetra Headers */
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"
#include "AztecOO.h"


#ifdef HAVE_ML_MATLAB
/* Needed for MEX */
#include "mex.h"

using namespace Teuchos;
using namespace ML_Epetra;
using namespace MLAPI;
using namespace std;

extern void _main();

/* Macros */
#define MAX(x,y) ((x)>(y)?(x):(y))
#define FLOOR(x) ((int)(x))
#define ISINT(x) (((x-(int)(x))<1e-15*(x))?true:false)
#define IS_FALSE 0
#define IS_TRUE  1
#define MLMEX_ERROR -1

/* Mode Info */
typedef enum {MODE_SETUP=0, MODE_SOLVE, MODE_CLEANUP, MODE_STATUS, MODE_AGGREGATE} MODE_TYPE;

/* MLMEX Teuchos Parameters*/
#define MLMEX_INTERFACE "mlmex: interface"

/* Default values */
#define MLMEX_DEFAULT_LEVELS 10
#define MLMEX_DEFAULT_NUMPDES 1
#define MLMEX_DEFAULT_ADAPTIVEVECS 0
#define MLMEX_DEFAULT_USEDEFAULTNS true

/* Debugging */
//#define VERBOSE_OUTPUT


/**************************************************************/
/**************************************************************/
/**************************************************************/

#ifdef TESTING_FUNCTIONS_ONLY
/* Printout function for testing */
void csc_print(int n,int *rowind,int* colptr, double* vals){
  int i,j;
  for(i=0;i<n;i++)
    for(j=colptr[i];j<colptr[i+1];j++)
      mexPrintf("%d %d %20.16e\n",rowind[j],i,vals[j]);  
}/*end if*/
#endif


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* The big classes of data */
class ml_data_pack;
class ml_data_pack{
public:
  ml_data_pack();
  virtual ~ml_data_pack();

  /* setup - sets up an ml_data_pack object
     Parameters:
     N       - Number of unknowns [I]
     rowind  - Row indices of matrix (CSC format) [I]
     colptr  - Column indices of matrix (CSC format) [I]
     vals    - Nonzero values of matrix (CSC format) [I]
     Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  virtual int setup(int N,int* rowind,int* colptr, double* vals)=0;

  /* status - reports (to stdout) the status of the object
     Returns: IS_TRUE 
  */
  virtual int status()=0;

  /* solve - Given two Teuchos lists, one in the ml_data_pack, and one of
     solve-time options, this routine calls the relevant solver and returns the solution.
     Parameters:
     TPL     - Teuchos list of solve-time options [I]
     N       - Number of unknowns [I]
     b       - RHS vector [I]
     x       - solution vector [O]
     returns: IS_TRUE if sucessful, IS_FALSE otherwise.
  */
  virtual int solve(Teuchos::ParameterList* TPL, int N, double*b, double*x)=0;  
public:
  int id;
  Teuchos::ParameterList *List;
  double operator_complexity;  
  ml_data_pack *next;
};

class mlapi_data_pack:public ml_data_pack{
public:
  mlapi_data_pack();
  ~mlapi_data_pack();

  /* setup - sets up an ml_data_pack object
     Parameters:
     N       - Number of unknowns [I]
     rowind  - Row indices of matrix (CSC format) [I]
     colptr  - Column indices of matrix (CSC format) [I]
     vals    - Nonzero values of matrix (CSC format) [I]
     Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  int setup(int N,int* rowind,int* colptr, double* vals);

  /* status - reports (to stdout) the status of the object
     Returns: IS_TRUE 
  */
  int status();


  /* solve - Given two Teuchos lists, one in the ml_data_pack, and one of
     solve-time options, this routine calls the relevant solver and returns the solution.
     Parameters:
     TPL     - Teuchos list of solve-time options [I]
     N       - Number of unknowns [I]
     b       - RHS vector [I]
     x       - solution vector [O]
     returns: IS_TRUE if sucessful, IS_FALSE otherwise.
  */
  int solve(Teuchos::ParameterList* TPL, int N, double*b, double*x);  

private:
  Space * FineSpace;
  DistributedMatrix * A;
  MultiLevelAdaptiveSA *Prec;
};


class ml_epetra_data_pack:public ml_data_pack{
public:
  ml_epetra_data_pack();
  ~ml_epetra_data_pack();

  /* setup - sets up an ml_data_pack object
     Parameters:
     N       - Number of unknowns [I]
     rowind  - Row indices of matrix (CSC format) [I]
     colptr  - Column indices of matrix (CSC format) [I]
     vals    - Nonzero values of matrix (CSC format) [I]
     Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
  */
  int setup(int N,int* rowind,int* colptr, double* vals);

  /* status - reports (to stdout) the status of the object
     Returns: IS_TRUE 
  */
  int status();

  /* solve - Given two Teuchos lists, one in the ml_data_pack, and one of
     solve-time options, this routine calls the relevant solver and returns the solution.
     Parameters:
     TPL     - Teuchos list of solve-time options [I]
     N       - Number of unknowns [I]
     b       - RHS vector [I]
     x       - solution vector [O]
     returns: IS_TRUE if sucessful, IS_FALSE otherwise.
  */
  int solve(Teuchos::ParameterList* TPL, int N, double*b, double*x);  

private:
  Epetra_Comm *Comm;
  Epetra_Map *Map;
  Epetra_CrsMatrix * A;
  MultiLevelPreconditioner *Prec;
};


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* ML data pack list */
class ml_data_pack_list;
class ml_data_pack_list{
public:
  ml_data_pack_list();
  ~ml_data_pack_list();

  /* add - Adds an ML_DATA_PACK to the list.
     Parameters:
     D       - The ML_DATA_PACK. [I]
     Returns: problem id number of D
  */
  int add(ml_data_pack *D);

  /* find - Finds problem by id
     Parameters:
     id      - ID number [I]
     Returns: pointer to ML_DATA_PACK matching 'id', if found, NULL if not
     found.
  */
  ml_data_pack* find(int id);

  /* remove - Removes problem by id
     Parameters:
     id      - ID number [I]
     Returns: IS_TRUE if remove was succesful, IS_FALSE otherwise
  */
  int remove(int id);

  /* size - Number of stored problems
     Returns: num_probs
  */
  int size();  

  /* Returns the status of all members of the list
     Returns IS_TRUE
  */  
  int status_all();

  
protected:
  int num_probs;
  /* Note: This list is sorted */
  ml_data_pack *L;
};


/**************************************************************/
/**************************************************************/
/**************** ml_data_pack class functions ****************/
/**************************************************************/
/**************************************************************/
ml_data_pack::ml_data_pack():id(MLMEX_ERROR),List(NULL),operator_complexity(MLMEX_ERROR),next(NULL){}
ml_data_pack::~ml_data_pack(){if(List) delete List;}


/**************************************************************/
/**************************************************************/
/*************** mlapi_data_pack class functions **************/
/**************************************************************/
/**************************************************************/
mlapi_data_pack::mlapi_data_pack():ml_data_pack(),FineSpace(NULL),A(NULL),Prec(NULL){}

mlapi_data_pack::~mlapi_data_pack(){
   if(FineSpace) delete FineSpace;
   if(A)         delete A;
   if(Prec)      delete Prec;
}/*end destructor*/


/* mlapi_data_pack::status - This function does a status query on the
   MLAPI_DATA_PACK passed in.
   Returns: IS_TRUE
*/
int mlapi_data_pack::status(){
  mexPrintf("**** Problem ID %d [MLAPI] ****\n",id);
  if(A) mexPrintf("Matrix: %dx%d w/ %d nnz\n",A->NumGlobalRows(),A->NumGlobalCols(),A->NumGlobalNonzeros()); 
  mexPrintf(" Operator complexity = %e\n",operator_complexity);
  if(List){mexPrintf("Parameter List:\n");List->print(cout,1);}
  mexPrintf("\n");
  return IS_TRUE;
}/*end status*/
                                      
/**************************************************************/
/**************************************************************/
/**************************************************************/
/* mlapi_data_pack::setup - This function does the setup phase for MLAPI, pulling
   key parameters from the Teuchos list, and calling the aggregation routines
   Parameters:
   N       - Number of unknowns [I]
   rowind  - Row indices of matrix (CSC format) [I]
   colptr  - Column indices of matrix (CSC format) [I]
   vals    - Nonzero values of matrix (CSC format) [I]
   Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
*/
int mlapi_data_pack::setup(int N,int* rowind,int* colptr, double* vals){

  int i,j;
  /* Initialize the workspace and set the output level */
  Init();  

  /* Define the space for fine level vectors and operators */
  FineSpace=new Space(N);
  
  /* Do the matrix assembly --- Someone should really come up with a better way
     of doing this.  Argh! */
  A=new DistributedMatrix(*FineSpace,*FineSpace); 
  for(i=0;i<N;i++)
    for(j=colptr[i];j<colptr[i+1];j++)
      (*A)(rowind[j],i)=vals[j];
  A->FillComplete();

  /* Pull key options from Teuchos */
  int numpdes=List->get("PDE equations",MLMEX_DEFAULT_NUMPDES);
  int cutoff=List->get("max levels",MLMEX_DEFAULT_LEVELS);
  int adaptivevecs=List->get("additional candidates",MLMEX_DEFAULT_ADAPTIVEVECS);
  bool UseDefaultNullSpace=List->get("use default null space", MLMEX_DEFAULT_USEDEFAULTNS);

  /* Allocate Memory */
  Prec=new MultiLevelAdaptiveSA(*A,*List,numpdes,cutoff);
    
  /* Build the Heirarchy */
  if(adaptivevecs>0) Prec->AdaptCompute(UseDefaultNullSpace,adaptivevecs);    
  else Prec->Compute();

  operator_complexity=Prec->GetComplexity();
  printf("Smoothed Aggregation: operator complexity = %e\n",operator_complexity);
  
  Finalize();
  return IS_TRUE;
}/*end setup*/


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* mlapi_data_pack::solve - Given two Teuchos lists, one in the MLAPI_DATA_PACK, and one of
   solve-time options, this routine calls the relevant solver and returns the solution.
   Parameters:
   TPL     - Teuchos list of solve-time options [I]
   N       - Number of unknowns [I]
   b       - RHS vector [I]
   x       - solution vector [O]
   Returns: IS_TRUE if solve was succesful, IS_FALSE otherwise
*/
int mlapi_data_pack::solve(Teuchos::ParameterList *TPL, int N, double*b, double*x){
  int i;
  Init();   
  MultiVector LHS(A->GetDomainSpace());
  MultiVector RHS(A->GetRangeSpace());

  /* Create the new Teuchos List */
  Teuchos::ParameterList tmp=*List;
  tmp.setParameters(*TPL);
  
  /* Fill LHS/RHS */
  LHS=0;
  for(i=0;i<N;i++) RHS(i)=b[i];

#ifdef VERBOSE_OUTPUT
  tmp.print(cout,0,true);
#endif
  
  /* Do the solve */
  Krylov(*A,LHS,RHS,*Prec,tmp);

  /* Fill Solution */
  for(i=0;i<N;i++) x[i]=LHS(i);
  
  Finalize();
  return IS_TRUE;
}/*end solve*/


/**************************************************************/
/**************************************************************/
/************* ml_epetra_data_pack class functions ************/
/**************************************************************/
/**************************************************************/
ml_epetra_data_pack::ml_epetra_data_pack():ml_data_pack(),Comm(NULL),Map(NULL),A(NULL),Prec(NULL){}
ml_epetra_data_pack::~ml_epetra_data_pack(){
  if(Comm) delete Comm;
  if(Map)  delete Map;
  if(A)    delete A;
  if(Prec) delete Prec;
}/*end destructor*/

/* ml_epetra_data_pack_status - This function does a status query on the
   ML_EPETRA_DATA_PACK passed in.
   Returns: IS_TRUE
*/
int ml_epetra_data_pack::status(){
  mexPrintf("**** Problem ID %d [ML_Epetra] ****\n",id);
  if(A) mexPrintf("Matrix: %dx%d w/ %d nnz\n",A->NumGlobalRows(),A->NumGlobalCols(),A->NumMyNonzeros()); 
  mexPrintf(" Operator complexity = %e\n",operator_complexity);
  if(List){mexPrintf("Parameter List:\n");List->print(cout,1);}
  mexPrintf("\n");
  return IS_TRUE;
}/*end status*/


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* ml_epetra_data_pack::setup - This function does the setup phase for ML_Epetra, pulling
   key parameters from the Teuchos list, and calling the aggregation routines
   Parameters:
   N       - Number of unknowns [I]
   rowind  - Row indices of matrix (CSC format) [I]
   colptr  - Column indices of matrix (CSC format) [I]
   vals    - Nonzero values of matrix (CSC format) [I]
   Returns: IS_TRUE if setup was succesful, IS_FALSE otherwise
*/
int ml_epetra_data_pack::setup(int N,int* rowind,int* colptr, double* vals){
  int i,j;
  int *rnz;
  
  /* Nonzero counts for Epetra */
  rnz=new int[N];
  for(i=0;i<N;i++) rnz[i]=rowind[i+1]-rowind[i];  
  
  /* Epetra Setup */
  Comm= new Epetra_SerialComm;
  Map=new Epetra_Map(N,0,*Comm);
  A=new Epetra_CrsMatrix(Copy,*Map,rnz);
  
  /* Do the matrix assembly */
  for(i=0;i<N;i++)
    for(j=colptr[i];j<colptr[i+1];j++)
      A->InsertGlobalValues(rowind[j],1,&vals[j],&i);
  //NTS: Redo with block assembly, remembering to transpose
  A->FillComplete();

  /* Allocate Memory */
  Prec=new MultiLevelPreconditioner(*A, *List,false);  
  
  /* Build the Heirarchy */
  Prec->ComputePreconditioner();

  /* Pull Operator Complexity */
  operator_complexity = Prec->GetML_Aggregate()->operator_complexity / Prec->GetML_Aggregate()->fine_complexity;
       
  return IS_TRUE;
}/*end setup*/


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* ml_epetra_data_pack::solve - Given two Teuchos lists, one in the ml_epetra_data_pack, and one of
   solve-time options, this routine calls the relevant solver and returns the solution.
   Parameters:
   TPL     - Teuchos list of solve-time options [I]
   N       - Number of unknowns [I]
   b       - RHS vector [I]
   x       - solution vector [O]
   Returns: IS_TRUE if solve was succesful, IS_FALSE otherwise
*/
int ml_epetra_data_pack::solve(Teuchos::ParameterList *TPL, int N, double*b, double*x){
  int i;
  Epetra_Vector LHS(*Map);
  Epetra_Vector RHS(*Map);

  /* Fill LHS/RHS */
  LHS.PutScalar(0);
  for(i=0;i<N;i++) RHS[i]=b[i];

  /* Create the new Teuchos List */
  Teuchos::ParameterList tmp=*(List);
  tmp.setParameters(*TPL);
#ifdef VERBOSE_OUTPUT
  tmp.print(cout,0,true);
#endif
  
  /* Define Problem / Solver */
  Epetra_LinearProblem Problem(A, &LHS, &RHS);
  AztecOO solver(Problem);
  solver.SetPrecOperator(Prec);

  /* Get solver options from Teuchos list */
  int    NumIters = tmp.get("krylov: max iterations", 1550);
  double Tol      = tmp.get("krylov: tolerance", 1e-9);
  string type     = tmp.get("krylov: type", "gmres");
  int    output   = tmp.get("krylov: output level",10);
  string conv     = tmp.get("krylov: conv", "r0");  

  /* Set solver options - Solver type*/
  if (type == "cg") solver.SetAztecOption(AZ_solver, AZ_cg);
  else if (type == "cg_condnum") solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  else if (type == "gmres") solver.SetAztecOption(AZ_solver, AZ_gmres);
  else if (type == "gmres_condnum") solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
  else if (type == "fixed point") solver.SetAztecOption(AZ_solver, AZ_fixed_pt);
  //  else ML_THROW("krylov: type has incorrect value (" + type + ")", -1);
  
  /* Set solver options - Convergence Criterion*/
  if(conv == "r0") solver.SetAztecOption(AZ_conv,AZ_r0);
  else if(conv == "rhs") solver.SetAztecOption(AZ_conv,AZ_rhs);
  else if(conv == "Anorm") solver.SetAztecOption(AZ_conv,AZ_Anorm);
  else if(conv == "noscaled") solver.SetAztecOption(AZ_conv,AZ_noscaled);
  else if(conv == "sol") solver.SetAztecOption(AZ_conv,AZ_sol);
  //  else ML_THROW("krylov: conv has incorrect value (" + conv + ")",-1);

  /* Set solver options - other */
  solver.SetAztecOption(AZ_output, output);

  /* Do the solve */
  solver.Iterate(NumIters, Tol);

  /* Fill Solution */
  for(i=0;i<N;i++) x[i]=LHS[i];

  return IS_TRUE;
}/*end solve*/



/**************************************************************/
/**************************************************************/
/***********  ml_data_pack list class functions ***************/
/**************************************************************/
/**************************************************************/

ml_data_pack_list::ml_data_pack_list():num_probs(0),L(NULL){}
ml_data_pack_list::~ml_data_pack_list(){
  ml_data_pack *old;
  while(L!=NULL){
    old=L; L=L->next;
    delete old;
    /*NTS: Let's hope inheritance works and the correct destructor is called.  I
      have my doubts.*/
  }/*end while*/
}/*end destructor*/


/* add - Adds an ml_data_pack to the list.
   Parameters:
   D       - The ml_data_pack. [I]
   Returns: problem id number of D
*/
int ml_data_pack_list::add(ml_data_pack *D){
  int idx_prev=0;
  ml_data_pack *iprev, *icurr;
  
  if(!L){L=D;D->next=NULL;D->id=1;}
  else {
    /* Find the first numbering gap + add in D appropriately */
    for(iprev=NULL,icurr=L; icurr && icurr->id-idx_prev==1; iprev=icurr,idx_prev=iprev->id,icurr=icurr->next);
    if(!iprev) L=D;
    else iprev->next=D;
    D->id=idx_prev+1;
    D->next=icurr;
    num_probs++;
  }/*end if*/  
  return D->id;
}/*end add*/


/* find - Finds problem by id
   Parameters:
   id      - ID number [I]
   Returns: pointer to ml_data_pack matching 'id', if found, NULL if not
   found.
*/
ml_data_pack* ml_data_pack_list::find(int id){
  ml_data_pack *rv;
  for(rv=L;rv && rv->id<id;rv=rv->next);
  if(rv && rv->id==id) return rv;
  else return NULL;  
}/*end find*/


/* remove - Removes problem by id
   Parameters:
   id      - ID number [I]
   Returns: IS_TRUE if remove was succesful, IS_FALSE otherwise
*/
int ml_data_pack_list::remove(int id){
  ml_data_pack *iprev, *icurr;
  if(!L) return IS_FALSE;
  for(iprev=NULL,icurr=L; icurr && icurr->id<id; iprev=icurr,icurr=icurr->next);

  if(!icurr || icurr->id!=id) return IS_FALSE;
  else{
    if(!iprev) L=icurr->next;
    else iprev->next=icurr->next;
    delete icurr;
    /*NTS: Hope that inheritance does its magic and that the right destructor is
      called */
    return IS_TRUE;
  }/*end else*/    
}/*end remove*/


/* size - Number of stored problems
   Returns: num_probs
*/
int ml_data_pack_list::size(){return num_probs;}


/* Returns the status of all members of the list
   Returns IS_TRUE
*/  
int ml_data_pack_list::status_all(){
  ml_data_pack *curr;
  for(curr=L;curr;curr=curr->next)
    curr->status();
  return IS_TRUE;
}/*end status_all */



/**************************************************************/
/**************************************************************/
/********************* Utility Functions **********************/
/**************************************************************/
/**************************************************************/
/* sanity_check - sanity checks the first couple of arguements and returns the
   program mode.
   Parameters:
   nrhs    - Number of program inputs [I]
   prhs    - The problem inputs [I]
   Return value: Which mode to run the program in.  
*/

MODE_TYPE sanity_check(int nrhs, const mxArray *prhs[]){
  MODE_TYPE rv;
  double *modes;
  /* Check for mode */
  if(nrhs==0)
    mexErrMsgTxt("Error: Invalid Inputs\n");
  
  /* Pull mode data from 1st Input */
  modes=mxGetPr(prhs[0]);

  switch ((MODE_TYPE)modes[0]){
  case MODE_SETUP:
    if(nrhs>1&&mxIsSparse(prhs[1])) rv=MODE_SETUP;
    else mexErrMsgTxt("Error: Invalid input for setup\n");    
    break;
  case MODE_SOLVE:
    if(nrhs>3&&mxIsNumeric(prhs[1])&&mxIsSparse(prhs[2])&&mxIsNumeric(prhs[3])) rv=MODE_SOLVE;
    else mexErrMsgTxt("Error: Invalid input for solve\n");
    break;
  case MODE_CLEANUP:
    if(nrhs==1 || nrhs==2) rv=MODE_CLEANUP;
    else mexErrMsgTxt("Error: Extraneous args for cleanup\n");
    break;
  case MODE_STATUS:
    if(nrhs==1 || nrhs==2) rv=MODE_STATUS;
    else mexErrMsgTxt("Error: Extraneous args for status\n");
    break;    
  default:
    mexErrMsgTxt("Error: Invalid input mode\n");
  };
  return rv;
}/*end sanity_check*/


/**************************************************************/
/**************************************************************/
/**************************************************************/
/* build_teuchos_list - takes the inputs (barring the solver mode and
  matrix/rhs) and turns them into a Teuchos list for use by MLAPI.
   Parameters:
   nrhs    - Number of program inputs [I]
   prhs    - The problem inputs [I]
   Return value: Teuchos list containing all parameters passed in by the user.
*/
Teuchos::ParameterList* build_teuchos_list(int nrhs,const mxArray *prhs[]){
  Teuchos::ParameterList* TPL=new Teuchos::ParameterList;
  mxClassID cid;
  int i,M,N, *opt_int;
  char *option_name,*opt_char;
  double *opt_float;
  string opt_str;  
  for(i=0;i<nrhs;i+=2){
    if(i==nrhs-1 || !mxIsChar(prhs[i]))
      mexErrMsgTxt("Error: Input options are not in ['parameter',value] format!\n");

    /* What option are we setting? */
    option_name=mxArrayToString(prhs[i]);

    /* Pull relevant info the the option value */
    cid=mxGetClassID(prhs[i+1]);
    M=mxGetM(prhs[i+1]);
    N=mxGetN(prhs[i+1]);
#ifdef VERBOSE_OUTPUT
    mexPrintf("[%d] M=%d N=%d\n",i,M,N);
#endif

    /* Add to the Teuchos list */
    switch(cid){
    case mxCHAR_CLASS:
      opt_char=mxArrayToString(prhs[i+1]);
      opt_str=opt_char;
#ifdef VERBOSE_OUTPUT
      mexPrintf("[%s] String Found: %s\n",option_name,opt_char);
#endif
      TPL->set(option_name, opt_str);
      mxFree(opt_char);
      break;
    case mxDOUBLE_CLASS:
    case mxSINGLE_CLASS:
      //NTS: Does not deal with complex args
      opt_float=mxGetPr(prhs[i+1]);
      if(M==1 && N==1 && ISINT(opt_float[0])) {     
#ifdef VERBOSE_OUTPUT
        mexPrintf("[%s] Float(Int) Found!\n",option_name);
#endif
        TPL->set(option_name, (int)opt_float[0]);
      }/*end if*/
      else if(M==1 && N==1){
#ifdef VERBOSE_OUTPUT
        mexPrintf("[%s] Float Found!\n",option_name);
#endif
        TPL->set(option_name, opt_float[0]);
      }/*end if*/
      else{
#ifdef VERBOSE_OUTPUT
        mexPrintf("[%s] Float Found!\n",option_name);
#endif
        TPL->set(option_name, opt_float);
      }/*end else*/
      break;
    case mxLOGICAL_CLASS:
#ifdef VERBOSE_OUTPUT      
      mexPrintf("[%s] Logical Found!\n",option_name);
#endif
      if(M==1 && N==1) TPL->set(option_name, mxIsLogicalScalarTrue(prhs[i+1]));
      else TPL->set(option_name,mxGetLogicals(prhs[i+1]));
      //NTS: The else probably doesn't work.
      break;
    case mxINT8_CLASS:
    case mxUINT8_CLASS:
    case mxINT16_CLASS:
    case mxUINT16_CLASS:
    case mxINT32_CLASS:
    case mxUINT32_CLASS:

#ifdef VERBOSE_OUTPUT
      mexPrintf("[%s] Int Found!\n",option_name);
#endif
      opt_int=(int*)mxGetData(prhs[i+1]);
      if(M==1 && N==1) TPL->set(option_name, opt_int[0]);      
      else TPL->set(option_name, opt_int);      
      break;
      // NTS: 64-bit ints will break on a 32-bit machine.  We
      // should probably detect machine type, or somthing, but that would
      // involve a non-trivial quantity of autoconf kung fu.
    case mxINT64_CLASS:
    case mxUINT64_CLASS:      
    case mxFUNCTION_CLASS:
    case mxUNKNOWN_CLASS:
    case mxCELL_CLASS:
    case mxSTRUCT_CLASS:      
    default:
      mexPrintf("Error parsing input option #%d: %s [type=%d]\n",i,option_name,cid);
      mexErrMsgTxt("Error: An input option is invalid!\n");      
    };        

    /* Free memory */
    mxFree(option_name);   
  }/*end for*/  

  /* Return */
  return TPL;
}/*end build_teuchos_list*/



/**************************************************************/
/**************************************************************/
/**************************************************************/
/* Calling syntax is done correctly by ml.m, please use that.
   The first rhs argument is always the program mode.  The other args are the
   ml.m parameter list in order (barring the mode arg).
*/

/* mexFunction is the gateway routine for the MEX-file. */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
  int i,*id,sz,rv,*rowind,*colptr;
  double *b, *x, *vals,*opts_array;
  int nr, nc,no;
  bool UseDefaultNullSpace;
  double tol,agg_thresh;
  string intf;
  Teuchos::ParameterList* List;
  mxClassID outclass;
  MODE_TYPE mode;
  static ml_data_pack_list* PROBS=NULL;
  ml_data_pack *D=NULL;
  
  /* Sanity Check Input */
  mode=sanity_check(nrhs,prhs);
  
  switch(mode){
  case MODE_SETUP:
    if(!PROBS) PROBS=new ml_data_pack_list;

    /* Pull Problem Size */
    nr=mxGetM(prhs[1]);
    nc=mxGetN(prhs[1]);

    /* Teuchos List*/
    if(nrhs>2) List=build_teuchos_list(nrhs-2,&(prhs[2]));
    else List=new Teuchos::ParameterList;

    /* Pick the Interface Type */
    intf=List->get(MLMEX_INTERFACE,"epetra");
    if(intf=="mlapi") D=new mlapi_data_pack();
    else D=new ml_epetra_data_pack();
    D->List=List;    
    
    /* Pull Matrix - CSC Format */
    vals=mxGetPr(prhs[1]);
    rowind=mxGetIr(prhs[1]);
    colptr=mxGetJc(prhs[1]);

    /* Construct the Heirarchy */
    D->setup(nr,rowind,colptr,vals);

    /* Add this problem to the list */
    rv=PROBS->add(D);

    /* Set return value(s) */
    plhs[0]=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    id=(int*)mxGetData(plhs[0]);id[0]=rv;
    if(nlhs>1) plhs[1]=mxCreateDoubleScalar(D->operator_complexity);
    
    /* Lock so we can keep the memory for the heirarchy */
    mexLock();
    break;

  case MODE_SOLVE:
    /* Are there problems set up? */
    if(!PROBS) mexErrMsgTxt("Error: No problems set up, cannot solve.\n");    

    /* Get the Problem Handle */
    id=(int*)mxGetData(prhs[1]);
    D=PROBS->find(id[0]);
    if(!D) mexErrMsgTxt("Error: Problem handle not allocated.\n");    
        
    /* Pull Problem Size */
    nr=mxGetM(prhs[2]);
    nc=mxGetN(prhs[2]);
    
    /* Pull RHS */
    b=mxGetPr(prhs[3]);

    /* Teuchos List*/
    if(nrhs>4) List=build_teuchos_list(nrhs-4,&(prhs[4]));
    else List=new Teuchos::ParameterList;
    
    /* Allocate Solution Space */
    plhs[0]=mxCreateDoubleMatrix(nr,1,mxREAL);
    x=mxGetPr(plhs[0]);
    
    /* Sanity Check Matrix / RHS */
    if(nr != nc || nr != mxGetM(prhs[2]))
      mexErrMsgTxt("Error: Size Mismatch in Input\n");

    /* Run Solver */  
    D->solve(List,nr,b,x);

    /* Cleanup */
    delete List;
    break;
    
  case MODE_CLEANUP:
    if(PROBS && nrhs==1){
      /* Cleanup all problems */
      sz=PROBS->size();
      for(i=0;i<sz;i++) mexUnlock();
      delete PROBS;PROBS=NULL;
      rv=1;
    }/*end if*/
    else if(PROBS && nrhs==2){
      /* Cleanup one problem */
      id=(int*)mxGetData(prhs[1]);
      rv=PROBS->remove(id[0]);
      if(rv) mexUnlock();
    }/*end elseif*/
    else rv=0;
    
    /* Set return value */
    plhs[0]=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    id=(int*)mxGetData(plhs[0]);id[0]=rv;    
    break;

  case MODE_STATUS:
    if(PROBS && nrhs==1){
      /* Status check on all problems */
      rv=PROBS->status_all();      
    }/*end if*/
    else if(PROBS && nrhs==2){
      /* Status check one problem */
      id=(int*)mxGetData(prhs[1]);
      D=PROBS->find(id[0]);
      if(!D) mexErrMsgTxt("Error: Problem handle not allocated.\n");      
      rv=D->status();      
    }/*end elseif*/
    else rv=0;
    
    /* Set return value */
    plhs[0]=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    id=(int*)mxGetData(plhs[0]);id[0]=rv;    
    break;
    
  default:
    mexErrMsgTxt("Error: Generic\n");
  };  
}/*end mexFunction*/

#else
#error "Do not have MATLAB"
#endif
#else
#error "Do not have MLAPI"
#endif
