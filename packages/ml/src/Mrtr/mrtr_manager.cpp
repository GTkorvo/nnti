/*
#@HEADER
# ************************************************************************
#
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifdef TRILINOS_PACKAGE

#include "mrtr_manager.H"
#include "EpetraExt_MatrixMatrix.h"  // for adding matrices
#include <EpetraExt_Transpose_RowMatrix.h>
#include "Epetra_Time.h"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
MOERTEL::Manager::Manager(Epetra_Comm& comm, int outlevel) :
outlevel_(outlevel),
comm_(comm),
dimensiontype_(MOERTEL::Manager::manager_none),
inputmap_(null),
inputmatrix_(null),
constraintsmap_(null),
D_(null),
M_(null),
saddlemap_(null),
saddlematrix_(null)
{
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
MOERTEL::Manager::~Manager()
{
  interface_.clear();
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 06/05|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const MOERTEL::Manager& man)
{ 
  man.Print();
  return (os);
}

/*----------------------------------------------------------------------*
 |  print all data                                           mwgee 06/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::Print() const
{ 
  comm_.Barrier();

  if (Comm().MyPID()==0)
  cout << "\n========================= Mortar Manager =========================\n\n";
  fflush(stdout);
  
  comm_.Barrier();

  map<int,RefCountPtr<MOERTEL::Interface> >::const_iterator curr;
  for (curr=interface_.begin(); curr!=interface_.end(); ++curr)
  {
    RefCountPtr<MOERTEL::Interface> inter = curr->second;
    if (inter==null)
    {
      cout << "***ERR*** MOERTEL::Manager::Print:\n"
           << "***ERR*** found NULL entry in map of interfaces\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    cout << *inter;
  }
  comm_.Barrier();
  if (inputmap_ != null)
  {
    if (comm_.MyPID() == 0)
    cout << "\n------------------ Input RowMap of original Problem ---------------\n";
    comm_.Barrier();
    cout << *inputmap_;
  }
  comm_.Barrier();
  if (constraintsmap_ != null)
  {
    if (comm_.MyPID() == 0)
    cout << "\n------------------ RowMap of Constraints ---------------\n";
    comm_.Barrier();
    cout << *constraintsmap_;
  }
  comm_.Barrier();
  if (D_ != null)
  {
    if (comm_.MyPID() == 0)
    cout << "\n------------------ Coupling Matrix D ---------------\n";
    comm_.Barrier();
    cout << *D_;
  }
  comm_.Barrier();
  if (M_ != null)
  {
    if (comm_.MyPID() == 0)
    cout << "\n------------------ Coupling Matrix M ---------------\n";
    comm_.Barrier();
    cout << *M_;
  }
  comm_.Barrier();

  if (Comm().MyPID()==0)
  cout << "\n========================= End Mortar Manager =========================\n\n";
  fflush(stdout);

  return true;
}

/*----------------------------------------------------------------------*
 |  Add an interface (public)                                mwgee 06/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::AddInterface(MOERTEL::Interface& interface) 
{
  if (!interface.IsComplete())
  {
    cout << "***ERR*** MOERTEL::Manager::AddInterface:\n"
         << "***ERR*** Cannot add segment as Complete() was not called before\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  RefCountPtr<MOERTEL::Interface> tmp = rcp(new MOERTEL::Interface(interface));
  interface_.insert(pair<int,RefCountPtr<MOERTEL::Interface> >(tmp->Id(),tmp));
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Store the rowmap of the underlying matrix from the application 07/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::SetInputMap(Epetra_Map* map)
{
  inputmap_ = rcp(new Epetra_Map(*map));
  return true;
}

/*----------------------------------------------------------------------*
 |  Set ptr to the uncoupled matrix on input                       07/05|
 |  Note that MOERTEL::Manager does not make a deep copy here by default   |
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::SetInputMatrix(Epetra_CrsMatrix* inputmatrix, bool DeepCopy)
{
  if (DeepCopy)
  {
    inputmatrix_ = rcp(new Epetra_CrsMatrix(*inputmatrix));
    return true;
  }
  else
  {
    inputmatrix_ = rcp(inputmatrix);
    inputmatrix_.release();
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------*
 |                                                                 08/05|
 |  modified version of the epetraext matrixmatrixadd                   |
 |  NOTE:                                                               |
 |  - A has to be FillComplete, B must NOT be FillComplete()            |
 *----------------------------------------------------------------------*/
int MOERTEL::Manager::MatrixMatrixAdd(const Epetra_CrsMatrix& A, bool transposeA,double scalarA,
                      Epetra_CrsMatrix& B,double scalarB )
{
  //
  //This method forms the matrix-matrix sum B = scalarA * op(A) + scalarB * B, where

  if (!A.Filled())
  {
     cout << "***ERR*** MOERTEL::Manager::MatrixMatrixAdd:\n"
          << "***ERR*** FillComplete was not called on A\n"
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
     exit(EXIT_FAILURE);
  }
  if (B.Filled())
  {
     cout << "***ERR*** MOERTEL::Manager::MatrixMatrixAdd:\n"
          << "***ERR*** FillComplete was called on B\n"
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
     exit(EXIT_FAILURE);
  }
  
  //explicit tranpose A formed as necessary
  Epetra_CrsMatrix* Aprime = 0;
  EpetraExt::RowMatrix_Transpose* Atrans = 0;
  if( transposeA )
  {
    Atrans = new EpetraExt::RowMatrix_Transpose(false);
    Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));
  }
  else
    Aprime = const_cast<Epetra_CrsMatrix*>(&A);
    
  //Loop over B's rows and sum into
  int MaxNumEntries = EPETRA_MAX( Aprime->MaxNumEntries(), B.MaxNumEntries() );
  int NumEntries;
  vector<int> Indices(MaxNumEntries);
  vector<double> Values(MaxNumEntries);

  int NumMyRows = A.NumMyRows();
  int Row, err;

  if( scalarA )
  {
    for( int i = 0; i < NumMyRows; ++i )
    {
      Row = Aprime->GRID(i);
      EPETRA_CHK_ERR(Aprime->ExtractGlobalRowCopy(Row,MaxNumEntries,NumEntries,&Values[0],&Indices[0]));
      if( scalarA != 1.0 )
        for( int j = 0; j < NumEntries; ++j ) Values[j] *= scalarA;
      if( B.Filled() ) {//Sum In Values
        err = B.SumIntoGlobalValues( Row, NumEntries, &Values[0], &Indices[0] );
        assert( err == 0 );
      }
      else {
        err = B.InsertGlobalValues( Row, NumEntries, &Values[0], &Indices[0] );
        assert( err == 0 || err == 1 );
      }
    }
  }

  Indices.clear();
  Values.clear();

  if( Atrans ) delete Atrans;

  return(0);
}

/*----------------------------------------------------------------------*
 |  Choose dofs for lagrange multipliers (private)           mwgee 07/05|
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_Map> MOERTEL::Manager::LagrangeMultiplierDofs()
{
  if (inputmap_==null)
  {
    cout << "***ERR*** MOERTEL::Manager::LagrangeMultiplierDofs:\n"
         << "***ERR*** inputmap==NULL, Need to set an input-rowmap first\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return null;
  }
  
  // find the largest row number in inputmap
  int maxinputGID = inputmap_->MaxAllGID();

  // start numbering lagrange multipliers from maxinputGID+1
  int minLMGID = maxinputGID+1;
  int maxLMGID = 0;
  
  int length = 0;
  
  // 2D problem:
  // loop interfaces and set LM dofs on the slave side for all slave nodes
  // that have a projection
  // 3D problem:
  // loop interfaces and set LM dofs on the slave side for all slave nodes
  // that have an integrated row of D and M projection

  // start with minLMGID and return maxLMGID+1 on a specific interface
  // Note this is collective for ALL procs
  map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    length -= minLMGID;
    maxLMGID = curr->second->SetLMDofs(minLMGID);
    if (!maxLMGID && maxLMGID!=minLMGID)
    {
      cout << "***ERR*** MOERTEL::Manager::LagrangeMultiplierDofs:\n"
           << "***ERR*** interface " << curr->second->Id() << " returned false\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return null;
    }
    minLMGID = maxLMGID;
    length += maxLMGID;
  }  

  vector<int> mylmids(length);
  int count=0;
  
  // get a vector of LM dof ids from each proc and each interface
  // and add it to the global one
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    RefCountPtr<MOERTEL::Interface> inter = curr->second;
    vector<int>* lmids = inter->MyLMIds();
    if (count+lmids->size() > mylmids.size())
      mylmids.resize(mylmids.size()+5*lmids->size());
    for (int i=0; i<(int)lmids->size(); ++i)
      mylmids[count++] = (*lmids)[i];
    delete lmids;
  }  
  mylmids.resize(count);
  int lsize = count;
  int gsize = 0;
  comm_.SumAll(&lsize,&gsize,1);
  
  // create the rowmap for the constraints
  // Note that this map contains the global communicator from the MOERTEL::Manager
  // NOT any interface local one
  RefCountPtr<Epetra_Map> map = rcp(new Epetra_Map(gsize,lsize,&(mylmids[0]),0,comm_));

  // tidy up
  mylmids.clear();
  
  return map;
}

/*----------------------------------------------------------------------*
 |  integrate all mortar interfaces                       (public) 11/05|
 |  Note: All processors have to go in here!                            |
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::Mortar_Integrate()
{
  //-------------------------------------------------------------------
  // check for problem dimension
  if (Dimension() == MOERTEL::Manager::manager_2D)
    return Mortar_Integrate_2D();
  else if (Dimension() == MOERTEL::Manager::manager_3D)
    return Mortar_Integrate_3D();
  else
  {
    cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
         << "***ERR*** problem dimension not set\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  return false;
}

/*----------------------------------------------------------------------*
 |  integrate all mortar interfaces in 2D                 (public) 07/05|
 |  Note: All processors have to go in here!                            |
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::Mortar_Integrate_2D()
{
  //-------------------------------------------------------------------
  // check for problem dimension
  if (Dimension() != MOERTEL::Manager::manager_2D)
  {
    cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
         << "***ERR*** problem dimension is not 2D?????\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  //-------------------------------------------------------------------
  // check whether we have interfaces
  int ninter = Ninterfaces();
  if (!ninter)
    return true;

  //-------------------------------------------------------------------
  // check whether we have an input map
  if (inputmap_==null)
  {
    cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
         << "***ERR*** inputmap==NULL, Need to set an input-rowmap first\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // check whether we have a mortar side chosen on each interface or 
  // whether we have to chose it here
  {
    map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
    bool foundit = true;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      int mside = curr->second->MortarSide();
      if (mside==-2)
      {
        foundit = false;
        break;
      }
    }
    if (!foundit) // we have to chose mortar sides ourself
      ChooseMortarSide();
  }  

  //-------------------------------------------------------------------
  // check whether functions have been set on interfaces
  // if not, check for functions flag and set them
  {
    bool foundit = true;
    map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      int nseg             = curr->second->GlobalNsegment();
      MOERTEL::Segment** segs = curr->second->GetSegmentView();
      for (int i=0; i<nseg; ++i)
        if (segs[i]->Nfunctions() < 1)
        {
          foundit = false;
          break;
        }
      delete [] segs;
      if (!foundit)
        curr->second->SetFunctionsFromFunctionTypes();
    }
  } 
  
  //-------------------------------------------------------------------
  // build projections for all interfaces
  {
    map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      bool ok  = curr->second->Project();
      if (!ok)
      {
        cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on projection\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }  

  //-------------------------------------------------------------------
  // this is probably the place to put detection of end segments
  // for each end segment, the order of the lagrange multiplier shape
  // function will be reduced by one
#if 0
  {
    map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      bool ok = curr->second->DetectEndSegmentsandReduceOrder();
      if (!ok)
      {
        cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false from DetectEndSegmentsandReduceOrder\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }
#endif

  //-------------------------------------------------------------------
  // choose dofs for lagrange multipliers and set them to slave nodes
  // build the rowmap for the coupling matrices M and D
  {
    constraintsmap_ = LagrangeMultiplierDofs();
    if (constraintsmap_==null)
    {
      cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
           << "***ERR*** LagrangeMultiplierDofs() returned NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }
    
  //-------------------------------------------------------------------
  // build the map for the saddle point problem
  {
    bool ok = BuildSaddleMap();
    if (!ok)
    {
      cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
           << "***ERR*** BuildSaddleMap() returned false\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }


  //-------------------------------------------------------------------
  // build the Epetra_CrsMatrix D and M
  D_ = rcp(new Epetra_CrsMatrix(Copy,*saddlemap_,5,false));
  M_ = rcp(new Epetra_CrsMatrix(Copy,*saddlemap_,40,false));

  //-------------------------------------------------------------------
  // integrate all interfaces
  {
    map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {  
      bool ok = curr->second->Mortar_Integrate(*D_,*M_);
      if (!ok)
      {
        cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on integration\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }
  
  //-------------------------------------------------------------------
  // call FillComplete() on M_ and D_ 
  D_->FillComplete(*saddlemap_,*saddlemap_);
  M_->FillComplete(*saddlemap_,*saddlemap_);
  
  //-------------------------------------------------------------------
  // print this
  if (OutLevel()>9)
    cout << *this;

  return true;
}

/*----------------------------------------------------------------------*
 |  integrate all mortar interfaces 3D case               (public) 11/05|
 |  Note: All processors have to go in here!                            |
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::Mortar_Integrate_3D()
{
  //-------------------------------------------------------------------
  // check for problem dimension
  if (Dimension() != MOERTEL::Manager::manager_3D)
  {
    cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
         << "***ERR*** problem dimension is not 3D?????\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  //-------------------------------------------------------------------
  // check whether we have interfaces
  int ninter = Ninterfaces();
  if (!ninter)
    return true;

  //-------------------------------------------------------------------
  // check whether we have an input map
  if (inputmap_==null)
  {
    cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
         << "***ERR*** inputmap==NULL, Need to set an input-rowmap first\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // check whether we have a mortar side chosen on each interface or 
  // whether we have to chose it here
  {
    map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
    bool foundit = true;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      int mside = curr->second->MortarSide();
      if (mside==-2)
      {
        foundit = false;
        break;
      }
    }
    if (!foundit) // we have to chose mortar sides ourself
      ChooseMortarSide();
  }  

  //-------------------------------------------------------------------
  // check whether functions have been set on interfaces
  // if not, check for functions flag and set them
  {
    bool foundit = true;
    map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      int nseg             = curr->second->GlobalNsegment();
      MOERTEL::Segment** segs = curr->second->GetSegmentView();
      for (int i=0; i<nseg; ++i)
        if (segs[i]->Nfunctions() < 1)
        {
          foundit = false;
          break;
        }
      delete [] segs;
      if (!foundit)
        curr->second->SetFunctionsFromFunctionTypes();
    }
  } 
  
  //-------------------------------------------------------------------
  // build normals for all interfaces
  {
    map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      bool ok = curr->second->BuildNormals();
      if (!ok)
      {
        cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on building normals\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }  

  //-------------------------------------------------------------------
  // this is probably the place to put detection of end segments
  // for each end segment, the order of the lagrange multiplier shape
  // function will be reduced by one
  // note that this is not yet implemented in 3D and is somehow
  // quite difficult to detect automatically. user input?
#if 1
  {
    map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      bool ok = curr->second->DetectEndSegmentsandReduceOrder();
      if (!ok)
      {
        cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false from DetectEndSegmentsandReduceOrder\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }
#endif

  //-------------------------------------------------------------------
  // integrate all interfaces
  // NOTE: 2D and 3D integration difer:
  // the 3D integration works without choosing lagrange multipliers
  // in advance. All integrated rows of D and M are stored as
  // scalars in the nodes
  {
    map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {  
      bool ok = curr->second->Mortar_Integrate();
      if (!ok)
      {
        cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on integration\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }

  //-------------------------------------------------------------------
  // choose dofs for lagrange multipliers and set them to slave nodes
  // build the rowmap for the coupling matrices M and D
  {
    constraintsmap_ = LagrangeMultiplierDofs();
    if (constraintsmap_==null)
    {
      cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
           << "***ERR*** LagrangeMultiplierDofs() returned NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }
  
  //-------------------------------------------------------------------
  // build the map for the saddle point problem
  {
    bool ok = BuildSaddleMap();
    if (!ok)
    {
      cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
           << "***ERR*** BuildSaddleMap() returned false\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }

  //-------------------------------------------------------------------
  // build the Epetra_CrsMatrix D and M
  D_ = rcp(new Epetra_CrsMatrix(Copy,*saddlemap_,50,false));
  M_ = rcp(new Epetra_CrsMatrix(Copy,*saddlemap_,40,false));

  //-------------------------------------------------------------------
  // now that we have all maps and dofs we can assemble from the nodes
  {
    map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {  
      bool ok = curr->second->Mortar_Assemble(*D_,*M_);
      if (!ok)
      {
        cout << "***ERR*** MOERTEL::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on assembly\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }

  //-------------------------------------------------------------------
  // call FillComplete() on M_ and D_ 
  D_->FillComplete(*saddlemap_,*saddlemap_);
  M_->FillComplete(*saddlemap_,*saddlemap_);
  
  //-------------------------------------------------------------------
  // print this
  if (OutLevel()>9)
    cout << *this;

  return true;
}


/*----------------------------------------------------------------------*
 | Create the saddle point problem map (private)             mwgee 11/05|
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::BuildSaddleMap()
{
  // check whether all interfaces are complete and integrated
  map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    if (curr->second->IsComplete() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::BuildSaddleMap:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not Complete()\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    if (curr->second->IsIntegrated() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::BuildSaddleMap:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not integrated yet\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }
  
  // check whether we have an inputmap
  if (inputmap_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::BuildSaddleMap:\n"
           << "***ERR*** No inputrowmap set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
  }
  
  // check whether we have a constraintsmap_
  if (constraintsmap_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::BuildSaddleMap:\n"
           << "***ERR*** onstraintsmap is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
  }
  
  // the saddle point problem rowmap is the inputmap + the constraintmap
  int numglobalelements = inputmap_->NumGlobalElements() +
                          constraintsmap_->NumGlobalElements();
  int nummyelements     = inputmap_->NumMyElements() +
                          constraintsmap_->NumMyElements();
  vector<int> myglobalelements(nummyelements);
  int count = 0;
  int* inputmyglobalelements = inputmap_->MyGlobalElements();
  for (int i=0; i<inputmap_->NumMyElements(); ++i)
    myglobalelements[count++] = inputmyglobalelements[i];
  int* constraintsmyglobalelements = constraintsmap_->MyGlobalElements();
  for (int i=0; i<constraintsmap_->NumMyElements(); ++i)
    myglobalelements[count++] = constraintsmyglobalelements[i];
  if (count != nummyelements)
  {
    cout << "***ERR*** MOERTEL::Manager::BuildSaddleMap:\n"
         << "***ERR*** Mismatch in dimensions\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  saddlemap_ = rcp(new Epetra_Map(numglobalelements,nummyelements,
                              &(myglobalelements[0]),0,comm_));
  myglobalelements.clear();
  
  return true;
}

/*----------------------------------------------------------------------*
 | Create the saddle point problem (public)                  mwgee 07/05|
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
Epetra_CrsMatrix* MOERTEL::Manager::MakeSaddleProblem()
{
  // check whether all interfaces are complete and integrated
  map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    if (curr->second->IsComplete() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not Complete()\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    if (curr->second->IsIntegrated() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not integrated yet\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
  }
  
  // check whether we have an inputmap
  if (inputmap_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** No inputrowmap set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
  }
  
  // check whether we have a constraintsmap_
  if (constraintsmap_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** onstraintsmap is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
  }
  
  // check for saddlemap_
  if (saddlemap_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** saddlemap_==NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
  }

  // check for inputmatrix
  if (inputmatrix_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** No inputmatrix set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
  }

  // check whether we have M and D matrices
  if (D_==null || M_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** Matrix M or D is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
  }
  
  // create a matrix for the saddle problem and fill it
  saddlematrix_ = rcp(new Epetra_CrsMatrix(Copy,*saddlemap_,90));

  // add values from inputmatrix
  MatrixMatrixAdd(*inputmatrix_,false,1.0,*saddlematrix_,0.0);
  
  // add values from D_
  MatrixMatrixAdd(*D_,false,1.0,*saddlematrix_,1.0);
  MatrixMatrixAdd(*D_,true,1.0,*saddlematrix_,1.0);  

  // add values from M_
  MatrixMatrixAdd(*M_,false,1.0,*saddlematrix_,1.0);
  MatrixMatrixAdd(*M_,true,1.0,*saddlematrix_,1.0);  

  saddlematrix_->FillComplete();

  return saddlematrix_.get();
}

/*----------------------------------------------------------------------*
 |                                                                 09/05|
 |  choose the mortar side                                              |
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::ChooseMortarSide()
{
  // find all interfaces
  vector<RefCountPtr<MOERTEL::Interface> > inter(Ninterfaces());
  map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
  curr=interface_.begin();
  int count = 0;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
      inter[count] = curr->second;
      ++count;
  }
  inter.resize(count);
  
  // call choice of the mortar side for all 1D interfaces
  bool ok = ChooseMortarSide(inter);

  // tidy up
  inter.clear();  

  return ok;
}


/*----------------------------------------------------------------------*
 |                                                                 09/05|
 |  choose the mortar side for 1D interfaces                            |
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::ChooseMortarSide(vector<RefCountPtr<MOERTEL::Interface> >& inter)
{
  // number of interfaces
  int ninter = inter.size();
  if (ninter < 2) 
  {
    if (inter[0]->MortarSide() != 0 && inter[0]->MortarSide() != 1)
      inter[0]->SetMortarSide(0);
    return true;
  }
  
  if (OutLevel()>5)
  {
    fflush(stdout);
    if (Comm().MyPID()==0)
    {
      cout << "---MOERTEL::Manager: start interface coloring:\n";
      for (int i=0; i<ninter; ++i)
        cout << "Interface " << inter[i]->Id() << " : Mortar Side : " << inter[i]->MortarSide() << endl;
      cout << "---MOERTEL::Manager: finding common nodes on interfaces:\n";
    }
    fflush(stdout);
  }
  Comm().Barrier();
  
  // time this process
  Epetra_Time time(Comm());
  time.ResetStartTime();
  
  
  // get a view of all nodes from all interfaces
  vector<MOERTEL::Node**> nodes(ninter);
  vector<int>          nnodes(ninter);
  for (int i=0; i<ninter; ++i)
  {
    nnodes[i] = inter[i]->GlobalNnode(); // returns 0 for procs not in lComm()
    nodes[i]  = inter[i]->GetNodeView(); // get vector of ALL nodes on inter i
  }
  
  // loop all interfaces
  for (int i=0; i<ninter; ++i)
  {
    // loop all nodes on that interface 
    // (procs not part of inter[i]->lComm() don't loop because they have nnodes[i]=0)
    for (int j=0; j<nnodes[i]; ++j)
    {
      // do nothing for node that does not belong to me
      if (inter[i]->NodePID(nodes[i][j]->Id()) != inter[i]->lComm()->MyPID())
        continue;
        
      // do nothing for a node that has been flagged cornernode before
      if (nodes[i][j]->IsCorner())
        continue;
        
      // the node we are currently looking for
      int actnodeid       = nodes[i][j]->Id();
      MOERTEL::Node* inode   = nodes[i][j];
      
      // search all other interfaces for this node
      for (int k=0; k<ninter; ++k)
      {
        if (k==i) 
          continue; // don't do the current interface
        
        MOERTEL::Node* knode = NULL;
        for (int l=0; l<nnodes[k]; ++l)
          if (actnodeid==nodes[k][l]->Id())
          {
            knode = nodes[k][l];
            break;
          }
        if (!knode) // node inode is not on interface k 
          continue;
          
        // found node actnodeid on interface i and k
        if (OutLevel()>5)
        {
          cout << "Node " << actnodeid << " on interfaces " << inter[i]->Id() << " and " << inter[k]->Id() << endl;
          fflush(stdout);
        }
        
        // flag that node on interfaces i and k as cornernode
        inode->SetCorner();
        knode->SetCorner();
      } // for (int k=0; k<ninter; ++k)
    } // for (int j=0; j< nnodes[i]; ++j)
  } // for (int i=0; i<ninter; ++i)

  
  // make the cornernode information redundant
    // loop all interfaces
  for (int i=0; i<ninter; ++i)
  {
    // if I'm not part of inter[i], continue
    if (!inter[i]->lComm())
      continue;
    for (int proc=0; proc<inter[i]->lComm()->NumProc(); ++proc)
    {
      // get the corner nodes I have found
      vector<int> bcast(4);
      int bsize = 0;
      if (proc==inter[i]->lComm()->MyPID())
        for (int j=0; j<nnodes[i]; ++j)
        {
          if (nodes[i][j]->IsCorner())
          {
            if ((int)bcast.size() <= bsize)
              bcast.resize(bsize+10);
            bcast[bsize] = j;
            ++bsize;  
          }
        } 
      inter[i]->lComm()->Broadcast(&bsize,1,proc);
      bcast.resize(bsize);
      inter[i]->lComm()->Broadcast(&(bcast[0]),bsize,proc);
      if (proc!=inter[i]->lComm()->MyPID())
        for (int k=0; k<bsize; ++k)
        {
          int j = bcast[k];
          nodes[i][j]->SetCorner();
        }
      bcast.clear();
    } // for (int proc=0; proc<Comm()->NumProc(); ++proc)
  } // for (i=0; i<ninter; ++i)

  
  // we have list of interfaces inter[i], now we need
  // to build list of cornernodes associated with each interface
  vector< vector<int> > cornernodes(ninter);
  for (int i=0; i<ninter; ++i) 
  {
    int bprocs = 0; int bprocr = 0;
    int bsize = 0;
    if (inter[i]->lComm())
      if (inter[i]->lComm()->MyPID()==0)
      {
        bprocs = Comm().MyPID();
        for (int j=0; j<nnodes[i]; ++j)
          if (nodes[i][j]->IsCorner())
          {
            if ((int)cornernodes[i].size() <= bsize)
              cornernodes[i].resize(bsize+5);
            cornernodes[i][bsize] = nodes[i][j]->Id();
            ++bsize;
          }
      }
    Comm().MaxAll(&bprocs,&bprocr,1);
    Comm().Broadcast(&bsize,1,bprocr);
    cornernodes[i].resize(bsize);
    Comm().Broadcast(&(cornernodes[i][0]),bsize,bprocr);
  }

  if (OutLevel()>5)
  {
    fflush(stdout);
    if (Comm().MyPID()==0)
      cout << "---MOERTEL::Manager: finding common nodes: " 
           << time.ElapsedTime() << " sec" << endl;
    fflush(stdout);
  }
  time.ResetStartTime();
  
  // we now color the interfaces in such a way that every corner node
  // has lagrange multipliers either never or once but never more then once
  // this exercice is a little tricky as each interface has it's own
  // intra-communicator where a bunch of procs but not necessarily all of them
  // belong to. Therefore some searching has to be done by procs that are 
  // member of a specific intra-comm, while some communication has to be
  // done using the global comm.
  // This loop probably does not scale well in parallel but I'm already happy
  // if it works at all!
  // The coloring allows for some interfaces to have a user-specified mortar side
  // It tries to make use of such a prescribed choice of side. But
  // when specifying a mortar side on one or more interfaces but not on all of them, 
  // there is no guarantee that the coloring works out even if the user chose
  // the mortar sides carefully so the coloring should/could work out!
  // This is due to the fact that there is some heuristics in here (e.g. line 926)
  // and there is more then one valid way to color the problem.
  vector<int> flags((2*ninter));
  vector<int> flagr((2*ninter));
  for (int i=0; i<ninter; ++i)
  {
     // check whether inter[i] has mortar side assigned
     // if not, go through all cornernodes and check dependency on
     // other interfaces
     // if there is a dependency, choose mortar side of inter[i] approp., 
     // if there is no dependency, choose something

     if (inter[i]->MortarSide() == -2)
     {
       // create a vector of flags as follows:
       // flags[0..ninter-1]        flags for side 0 of inter[i]
       // flags[ninter..2*ninter-1] flags for side 1 of inter[i]
       // flags[k]==  0 : no dependency
       // flags[k]==  1 : dependency but no side chosen yet on inter[k]
       // flags[k]==  2 : side should be chosen master/mortar side
       // flags[k]==  3 : side should be chosen slave/LM side
       for (int k=0; k<2*ninter; ++k) flags[k] = 0;
       
       // actnode[0] holds the active node id
       // actnode[1] holds the side this node is on on inter[i]  
       int actnodes[2]; 
       int actnoder[2]; 
       // loop cornernodes on inter[i]
       for (int j=0; j<(int)cornernodes[i].size(); ++j)
       {
         actnodes[0] = 0; actnodes[1] = 0;    
         if (inter[i]->lComm())
           if (inter[i]->lComm()->MyPID()==0)
           {
             actnodes[0] = cornernodes[i][j];
             actnodes[1] = inter[i]->GetSide(actnodes[0]);
           }
         Comm().MaxAll(actnodes,actnoder,2);
         
         // loop all other interfaces and check for actnoder[0]
         for (int k=0; k<ninter; ++k)
         {
           if (k==i) continue; // don't do inter[i]
           if (inter[k]->lComm())
             if (inter[k]->lComm()->MyPID()==0)
             {
               for (int l=0; l<(int)cornernodes[k].size(); ++l)
               {
                 if (actnoder[0]==cornernodes[k][l])
                 {
                   flags[k]        = 1; // found a dependency
                   flags[ninter+k] = 1; // found a dependancy
                   int mside = inter[k]->MortarSide();
                   if (mside != -2) // side has been chosen before on inter[k]
                   {
                     int nodeside = inter[k]->GetSide(actnoder[0]);
                     // found node on master/mortar side
                     // so set flags such that node is on slave on inter[i]
                     if (nodeside == mside) 
                     {
                       // actnode[1] has to become slave side on inter[i]
                       int sside_i = actnoder[1];
                       int mside_i = inter[k]->OtherSide(sside_i);
                       flags[sside_i*ninter+k] = 3;
                       flags[mside_i*ninter+k] = 2;
                     }
                     // found node on slave/LM side
                     // so set flags that node has to be on master on inter[i]
                     else if (nodeside == inter[k]->OtherSide(mside))
                     {
                       int mside_i = actnoder[1];
                       int sside_i = inter[k]->OtherSide(mside_i);
                       flags[sside_i*ninter+k] = 3;
                       flags[mside_i*ninter+k] = 2;
                     }
                     else
                     {
                       cout << "***ERR*** MOERTEL::Manager::ChooseMortarSide_2D:\n"
                            << "***ERR*** Matrix M or D is NULL\n"
                            << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
                       exit(EXIT_FAILURE);
                     }
                   }
                   break;
                 }
               } // for (int l=0; l<cornernodes[k].size(); ++l)
             }
         } // for (int k=0; k<ninter; ++k)
       } // for (int j=0; j<cornernodes[i].size(); ++j)

       // sum all flags
       Comm().MaxAll(&(flags[0]),&(flagr[0]),(2*ninter));

#if 0       
       cout << "Flags side 0:\n";
       for (int j=0; j<ninter; ++j) cout << "inter " << j << " flag " << flagr[j] << endl;
       cout << "Flags side 1:\n";
       for (int j=ninter; j<2*ninter; ++j) cout << "inter " << j << " flag " << flagr[j] << endl;
#endif
       
       // loop through flags and make a decision what to chose on inter[i]
       for (int k=0; k<ninter; ++k)
       {
         if (i==k) continue; 
         // I don't care what is chosen
         if (flagr[k]==0 && flagr[ninter+k]==0); 
         // i don't care either
         else if (flagr[k]==1 && flagr[ninter+k]==1);
         // side 0 has to be master, side 1 has to be slave
         else if (flagr[k]==2 && flagr[ninter+k]==3)
         {
           // check whether a side has been set before and this side is differen?
           if (inter[i]->MortarSide()==1)
           {
             cout << "***ERR*** MOERTEL::Manager::ChooseMortarSide_2D:\n"
                  << "***ERR*** want to set mortar side to 0 but has been set to 1 before\n"
                  << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
             exit(EXIT_FAILURE);
           }
           inter[i]->SetMortarSide(0);
         }
         // side 0 has to be slave, side 1 has to be master
         else if (flagr[k]==3 && flagr[ninter+k]==2)
         {
           // check whether a side has been set before and this side is differen?
           if (inter[i]->MortarSide()==0)
           {
             cout << "***ERR*** MOERTEL::Manager::ChooseMortarSide_2D:\n"
                  << "***ERR*** want to set mortar side to 1 but has been set to 0 before\n"
                  << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
             exit(EXIT_FAILURE);
           }
           inter[i]->SetMortarSide(1);
         }
         else
         {
           cout << "***ERR*** MOERTEL::Manager::ChooseMortarSide_2D:\n"
                << "***ERR*** unknown type of coloring flag: " << flags[k] << "\n"
                << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
           exit(EXIT_FAILURE);
         }
       }
       // Obviously we don't need to care what side inter[i] is going to have
       if (inter[i]->MortarSide()==-2)
         inter[i]->SetMortarSide(0);
         
       // now we loop through the flags again and color all interfaces,
       // that had a connection to inter[i] but did not care
       // (so they had flag[k]==1 and flag[ninter+k]==1)
       for (int k=0; k<ninter; ++k)
       {
         if (i==k) continue;
         if (flagr[k]==1 && flagr[ninter+k]==1)
         {
           if (inter[k]->MortarSide()!=-2)
           {
             cout << "***ERR*** MOERTEL::Manager::ChooseMortarSide_2D:\n"
                  << "***ERR*** weird, this interface is not supposed to already have a mortar side assigned\n"
                  << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
             exit(EXIT_FAILURE);
           }
           
           // loop through cornernodes of inter[i]
           for (int j=0; j<(int)cornernodes[i].size(); ++j)
           {
             actnodes[0] = 0; actnodes[1] = 0;    
             if (inter[i]->lComm())
               if (inter[i]->lComm()->MyPID()==0)
               {
                 actnodes[0] = cornernodes[i][j];
                 actnodes[1] = inter[i]->GetSide(actnodes[0]);
               }
             Comm().MaxAll(actnodes,actnoder,2);
             
             // intraproc 0 of inter[k] makes the decision on what to do
             int mortarside_s = -1;
             int mortarside_r;
             if (inter[k]->lComm())
               if (inter[k]->lComm()->MyPID()==0)
               {
                 // loop through inter[k]'s cornernode
                 for (int l=0; l<(int)cornernodes[k].size(); ++l)
                 {
                   if (actnoder[0]==cornernodes[k][l])
                   {
                     // get the mortar and nodal side of inter[i]
                     int mside_i   = inter[i]->MortarSide();
                     int nodeside_i = actnoder[1];
                     
                     // get the nodal side of inter[k]
                     int nodeside_k = inter[k]->GetSide(actnoder[0]);
                     
                     // the node is on the mortar side of i, 
                     // so put it on the slave side on k
                     if (mside_i == nodeside_i)
                     {
                       if (mortarside_s != -1)
                         if (mortarside_s != inter[k]->OtherSide(nodeside_k))
                         {
                           cout << "***ERR*** MOERTEL::Manager::ChooseMortarSide_2D:\n"
                                << "***ERR*** interface has a conflict that can not be resolved\n"
                                << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
                           exit(EXIT_FAILURE);
                         }
                       mortarside_s = inter[k]->OtherSide(nodeside_k);
                     }
                     // the node is on the slave side of i,
                     // so put it on the mortar side of k
                     else 
                     {
                       if (mortarside_s != -1)
                         if (mortarside_s != nodeside_k)
                         {
                           cout << "***ERR*** MOERTEL::Manager::ChooseMortarSide_2D:\n"
                                << "***ERR*** interface has a conflict that can not be resolved\n"
                                << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
                           exit(EXIT_FAILURE);
                         }
                       mortarside_s = nodeside_k;
                     }
                   }
                 } // for (l=0; l<cornernodes[k].size(); ++l)
               } // if (inter[k]->lComm()->MyPID()==0)
             Comm().MaxAll(&mortarside_s,&mortarside_r,1);
             if (mortarside_r != -1)
               inter[k]->SetMortarSide(mortarside_r);
           } // for (int j=0; j<cornernodes[i].size(); ++j)
           
         } // if (flagr[k]==1 && flagr[ninter+k]==1)
       } // for (k=0; k<ninter; ++k)

     } // if (inter[i]->MortarSide() == -2)
  } // for (int i=0; i<ninter; ++i)

  if (OutLevel()>5)
  {
    fflush(stdout);
    if (Comm().MyPID()==0)
      cout << "---MOERTEL::Manager: coloring interfaces : " 
           << time.ElapsedTime() << " sec" << endl;
    fflush(stdout);
  }

  if (OutLevel()>5)
  {
    fflush(stdout);
    if (Comm().MyPID()==0)
    {
      cout << "---MOERTEL::Manager: chosen mortar sides : \n";
      for (int i=0; i<ninter; ++i)
        cout << "Interface " << inter[i]->Id() << " : Mortar Side : " << inter[i]->MortarSide() << endl;
      cout << "---MOERTEL::Manager: end interface coloring\n";
    }
    fflush(stdout);
  }
  Comm().Barrier();

  // tidy up
  nnodes.clear();
  for (int i=0; i<ninter; ++i) 
    if (nodes[i]) delete [] nodes[i];
  nodes.clear();
  for (int i=0; i<ninter; ++i) 
    cornernodes[i].clear();
  cornernodes.clear();
  flags.clear();
  flagr.clear();

  return true;
}




















#endif // TRILINOS_PACKAGE
