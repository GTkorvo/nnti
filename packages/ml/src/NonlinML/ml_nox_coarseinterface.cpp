/*
#@HEADER
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
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
// ML-headers
#include "ml_common.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

// ----------   Includes   ----------
#include <ctime>
#include <iostream>
#include "ml_nox_coarseinterface.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            m.gee 12/04|
 *----------------------------------------------------------------------*/
ML_NOX::Nox_CoarseProblem_Interface::Nox_CoarseProblem_Interface(
                                    ML_NOX::Ml_Nox_Fineinterface& interface,
                                    int level,
                                    int plevel,
                                    Epetra_CrsMatrix** P,
                                    const Epetra_BlockMap* this_RowMap,
                                    int nlevel) 
: fineinterface_(interface)
{
  if (P==0)
  {
    cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::Nox_CoarseProblem_Interface:\n"
         << "**ERR**: Epetra_CrsMatrix** P is 0 in constructor on level " << level_ << "\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  int i;
  for (i=1; i<=level; i++)
  {
     if (P[i]==0)
     {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::Nox_CoarseProblem_Interface:\n"
             << "**ERR**: ptr to Prolongator level " << i << " is 0\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
     }
  }
  if (!this_RowMap)
  {
    cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::Nox_CoarseProblem_Interface:\n"
         << "**ERR**: Epetra_BlockMap* this_RowMap is 0 in constructor on level " << level_ << "\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }

  // store some data
  level_         = level;   // my level (ML_INCREASING only!!!)
  nlevel_        = nlevel;  // total number of levels
  ml_printlevel_ = plevel;  // set printing level
  nFcalls_       = 0;       // number of cals to the computeF function
  P_             = P;       // ptr to the vector of prolongators
  fbar_          = 0;       // the ptr to one FAS-vector
  fxbar_         = 0;       // the ptr to one FAS-vector
  isFASmodfied_  = false;
  this_RowMap_   = new Epetra_BlockMap(*this_RowMap);
  return;
}

/*----------------------------------------------------------------------*
 |  recreate the coarse interface (public)                   m.gee 12/04|
 |  returns false if recreation failed                                  |
 *----------------------------------------------------------------------*/
bool ML_NOX::Nox_CoarseProblem_Interface::recreate(int plevel,Epetra_CrsMatrix** P, 
                                                   const Epetra_BlockMap* this_RowMap) 
{
  if (P==0)
  {
    cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::Nox_CoarseProblem_Interface:\n";
    cout << "**ERR**: Epetra_CrsMatrix** P is 0 in constructor on level " << level_ << "\n";
    cout << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
  }
  int i;
  for (i=1; i<=level_; i++)
  {
     if (P[i]==0)
     {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::Nox_CoarseProblem_Interface:\n"
             << "**ERR**: ptr to Prolongator level " << i << " is 0\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
     }
  }
  // store some data
  ml_printlevel_ = plevel;  // set printing level
  nFcalls_       = 0;       // number of cals to the computeF function
  P_             = P;       // the new prolongators
  if (!this_RowMap)
  {
    cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::recreate:\n"
         << "**ERR**: Epetra_BlockMap* this_RowMap is 0 on level " << level_ << "\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  if (this_RowMap_) delete this_RowMap_;
  this_RowMap_   = new Epetra_BlockMap(*this_RowMap);
  return true;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            m.gee 12/04|
 *----------------------------------------------------------------------*/
ML_NOX::Nox_CoarseProblem_Interface::~Nox_CoarseProblem_Interface()
{ 
   // destroy the P-hierarchy, if it still exists
   destroyP();
   if (this_RowMap_) delete this_RowMap_;
   this_RowMap_ = 0;
   return; 
}

/*----------------------------------------------------------------------*
 |  restrict a vector from level current to                  m.gee 01/05|
 |  next coarser level next                                             |
 |  returns ptr to coarse Epetra_Vector                                 |
 |  NOTE: the calling routine is responsible for deleteing the          |
 |        returned Epetra_Vector                                        |
 *----------------------------------------------------------------------*/
Epetra_Vector* ML_NOX::Nox_CoarseProblem_Interface::restrict_to_next_coarser_level(
                                                                            Epetra_Vector* thisvec, 
                                                                            int current, int next)
{
   // check sanity of P_ - hierarchy
   if (P_==0)
   {
     cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_to_next_coarser_level:\n"
          << "**ERR**: P-hierarchy is NULL on level " << level_ << "\n"
          << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   int i;
   int max;
   if (current >= 1) max = current; else max = 1;
   for (i=max; i<=next; i++)
   {
      if (P_[i]==0)
      {
         cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_to_next_coarser_level:\n"
              << "**ERR**: ptr to Prolongator level " << i << " is NULL\n"
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
   }
   if (current+1 != next)
   {
     cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_to_next_coarser_level:\n"
          << "**ERR**: currently only impl. for consecutive levels\n"
          << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   
   // create a thisvec that matches the map of the restriction operator
   Epetra_Vector* xfineP = 0;
   if (current==0) // restrict from level 0 to level 1 is different from rest
   {
      // Note that the GIDs in xfine match those of the fineinterface and
      // might be different from those in P_[1]->RowMap.
      // The LIDs and the map match, so we have to copy xfine to xfineP
      // using LIDs
      xfineP = new Epetra_Vector(P_[1]->RowMap(),false);
      if (thisvec->MyLength() != xfineP->MyLength() || thisvec->GlobalLength() != xfineP->GlobalLength())
      {
          cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_to_next_coarser_level:\n"
               << "**ERR**: mismatch in dimension of thisvec and xfineP\n"
               << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      int mylength = thisvec->MyLength();
      for (i=0; i<mylength; i++)
         (*xfineP)[i] = (*thisvec)[i];      
   }
   else
   {
      xfineP = thisvec;
   }
   
   // Allocate matching coarse vector
   Epetra_Vector* cvec = new Epetra_Vector(P_[next]->OperatorDomainMap(),false);
   
   // restrict (transposed mutliply with prolongation)
   P_[next]->Multiply(true,*xfineP,*cvec);
   
   if (current==0)
      delete xfineP;
   xfineP=0;

   return cvec;
}

/*----------------------------------------------------------------------*
 |  prolong a vector from level current to                   m.gee 01/05|
 |  next coarser level next                                             |
 |  returns ptr to coarse Epetra_Vector                                 |
 |  NOTE: the calling routine is responsible for deleteing the          |
 |        returned Epetra_Vector                                        |
 *----------------------------------------------------------------------*/
Epetra_Vector* ML_NOX::Nox_CoarseProblem_Interface::prolong_to_this_level(Epetra_Vector* coarsevec, 
                                                                            int current, int next)
{
   // check sanity of P_ - hierarchy
   if (P_==0)
   {
     cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_to_next_coarser_level:\n"
          << "**ERR**: P-hierarchy is NULL on level " << level_ << "\n"
          << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   int i;
   int max = 0;
   if (current >= 1) max = current; else max = 1;
   for (i=max; i<=next; i++)
   {
      if (P_[i]==0)
      {
         cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_to_next_coarser_level:\n"
              << "**ERR**: ptr to Prolongator level " << i << " is NULL\n"
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
   }
   if (current+1 != next)
   {
     cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_to_next_coarser_level:\n"
          << "**ERR**: currently only impl. for consecutive levels\n"
          << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   
   // Allocate matching fine vector in ML-global indices
   Epetra_Vector* fvec = new Epetra_Vector(P_[next]->OperatorRangeMap(),false);
   
   // prolongate
   P_[next]->Multiply(false,*coarsevec,*fvec);
   
   // when prolongating to the finest level, global indexing there is different
   if (current==0)
   {
      Epetra_Vector* fine = new Epetra_Vector(fineinterface_.getMap(),false);
      if (fvec->MyLength() != fine->MyLength() || fvec->GlobalLength() != fine->GlobalLength())
      {
          cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_to_next_coarser_level:\n"
               << "**ERR**: mismatch in dimension of fvec and fine\n"
               << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      int mylength = fvec->MyLength();
      for (i=0; i<mylength; i++)
         (*fine)[i] = (*fvec)[i];
      delete fvec; fvec = 0;
      return fine;  
   }
   else
      return fvec;

   return NULL;
}

/*----------------------------------------------------------------------*
 |  restrict a fine level vector to this level               m.gee 12/04|
 |  NOTE: the calling routine is responsible for deleteing the          |
 |        returned Epetra_Vector                                        |
 *----------------------------------------------------------------------*/
Epetra_Vector* ML_NOX::Nox_CoarseProblem_Interface::restrict_fine_to_this(
                                               const Epetra_Vector& xfine)
{
   // on the finest level, just copy the vector
   if (level_==0)
   {
      Epetra_Vector* xcoarse = new Epetra_Vector(xfine);
      return xcoarse;
   }
   else // not on finest level
   {
      // check sanity of P_ - hierarchy
      if (P_==0)
      {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_fine_to_this:\n"
             << "**ERR**: P-hierarchy is NULL on level " << level_ << "\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      int i;
      for (i=1; i<=level_; i++)
      {
         if (P_[i]==0)
         {
            cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_fine_to_this:\n"
                 << "**ERR**: ptr to Prolongator level " << i << " is 0\n"
                 << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
         }
      }
      
      // Note that the GIDs in xfine match those of the fineinterface and
      // might be different from those in P_[1]->RowMap.
      // The LIDs and the map match, so we have to copy xfine to xfineP
      // using LIDs
      Epetra_Vector* xfineP = new Epetra_Vector(P_[1]->RowMap(),false);
      if (xfine.MyLength() != xfineP->MyLength() || xfine.GlobalLength() != xfineP->GlobalLength())
      {
          cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_fine_to_this:\n"
               << "**ERR**: mismatch in dimension of xfine and xfineP\n"
               << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
#if 0      
      cout << "xfine map\n";
      cout << xfine.Map();
      cout << "xfineP map\n";
      cout << xfineP->Map();
      exit(0);
#endif      
      int mylength = xfine.MyLength();
      for (i=0; i<mylength; i++)
         (*xfineP)[i] = xfine[i];      

      // loop from the finest level to this level and
      // apply series of restrictions (i.e transposed prolongations)
      Epetra_Vector* fvec = 0;
      for (i=0; i<level_; i++)
      {
         // allocate a vector matching level i+1
         Epetra_Vector* cvec = new Epetra_Vector(P_[i+1]->OperatorDomainMap(),false);
         
         // multiply
         if (i==0)
         {
            P_[i+1]->Multiply(true,*xfineP,*cvec);
            delete xfineP;
            xfineP = 0;
         }
         else
            P_[i+1]->Multiply(true,*fvec,*cvec);
         
         if (i>0)
            delete fvec;
         
         fvec = cvec;
      }
      if (xfineP) delete xfineP;
      xfineP = 0;
      return fvec;
   }
   return NULL;
}

/*----------------------------------------------------------------------*
 |  prolong a this coarse level vector to the finest level   m.gee 12/04|
 |  NOTE: the calling routine is responsible of deleteing the           |
 |        returned Epetra_Vector                                        |
 *----------------------------------------------------------------------*/
Epetra_Vector* ML_NOX::Nox_CoarseProblem_Interface::prolong_this_to_fine(
                                            const Epetra_Vector& xcoarse)
{
   // on the finest level, just copy the vector
   if (level_==0)
   {
      Epetra_Vector* xfine = new Epetra_Vector(xcoarse);
      return xfine;
   }
   else // we are not on the finest level
   {
      // check sanity of P_ - hierarchy
      if (P_==0)
      {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_fine_to_this:\n"
             << "**ERR**: P-hierarchy is NULL on level " << level_ << "\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      int i;
      for (i=1; i<=level_; i++)
      {
         if (P_[i]==0)
         {
            cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::restrict_fine_to_this:\n"
                 << "**ERR**: ptr to Prolongator level " << i << " is 0\n"
                 << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
         }
      }
      
      // loop from this coarsest level to the finest one
      // apply series of prolongations
      Epetra_Vector* cvec = 0;
      for (i=level_; i>0; i--)
      {
         // allocate a vector matching level i-1
         Epetra_Vector* fvec = new Epetra_Vector(P_[i]->OperatorRangeMap(),false);
         // multiply
         if (i==level_)
            P_[i]->Multiply(false,xcoarse,*fvec);
         else
            P_[i]->Multiply(false,*cvec,*fvec);
            
         if (i<level_)
            delete cvec;
            
         cvec = fvec;
      }
      // Note that the GIDs in cvec do NOT match those of the fineinterface as
      // they match the P_[1]->RangeMap.
      // The LIDs and the map match, so we have to copy cvec to xfine_fineinterface
      // using LIDs
      const Epetra_CrsGraph* graph_fineinterface = fineinterface_.getGraph();
      Epetra_Vector*   xfine_fineinterface = new Epetra_Vector(graph_fineinterface->RowMap(),false);
      if (cvec->MyLength() != xfine_fineinterface->MyLength() ||
          cvec->GlobalLength() != xfine_fineinterface->GlobalLength())
      {
          cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::prolong_this_to_fine:\n"
               << "**ERR**: mismatch in dimension of cvec and xfine_fineinterface\n"
               << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      int mylength = cvec->MyLength();
      for (i=0; i<mylength; i++)
         (*xfine_fineinterface)[i] = (*cvec)[i];
      delete cvec; cvec = 0;
      return xfine_fineinterface;
   }
   return NULL;
}

/*----------------------------------------------------------------------*
 |  evaluate nonlinear function (public, derived)            m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Nox_CoarseProblem_Interface::computeF(const Epetra_Vector& x, 
                                                   Epetra_Vector& F, 
                                                   const FillType fillFlag)
{
  bool err;
  int  ierr;
  ++nFcalls_;
  // check sanity of P_ - hierarchy
  if (P_==0)
  {
    cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeF:\n"
         << "**ERR**: P-hierarchy is NULL on level " << level_ << "\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  int i;
  for (i=1; i<=level_; i++)
  {
     if (P_[i]==0)
     {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeF:\n"
             << "**ERR**: ptr to Prolongator level " << i << " is 0\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
     }
  }

  // FIXME:: after intensive testing, this test might not be necessary
  Epetra_Vector* xcoarse  = 0;
#if 0
  bool samemap = this_RowMap_->PointSameAs(x.Map());
  if (samemap)
  {
#endif
     xcoarse = new Epetra_Vector(*this_RowMap_,false);
     xcoarse->Update(1.0,x,0.0);
#if 0
  }
  else
  {
     cout << "**WRN** Maps are not equal in\n"
          << "**WRN** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
     // this exporter exports from input map of x to the current graph_.RowMap
     Epetra_Export* exporter = new Epetra_Export(x.Map(),*this_RowMap_);     
     xcoarse = new Epetra_Vector(*this_RowMap_,false);
     ierr = xcoarse->Export(x,*exporter,Insert);
     if (ierr)
     {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeF:\n"
             << "**ERR**: export from x to xcoarse returned err=" << ierr <<"\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
     }
     delete exporter; exporter = 0;
  }
#endif  
  if (ml_printlevel_>9 && x.Comm().MyPID()==0)
  {
     cout << "ML (level " << level_ << "): Call no " << nFcalls_ << " to Nox_CoarseProblem_Interface::computeF\n";
     fflush(stdout);
  }     


  if (level_==0)
  {
     // create Ffine and xfine matching the fine interface
     const Epetra_CrsGraph* finegraph = fineinterface_.getGraph();
     Epetra_Vector*         Ffine     = new Epetra_Vector(finegraph->RowMap(),false);
     Epetra_Vector*         xfine     = new Epetra_Vector(finegraph->RowMap(),false);

     // FIXME:: after intensive testing, this test might not be necessary
#if 0
     samemap = xfine->Map().PointSameAs(xcoarse->Map());
     if (samemap)
     {
#endif
        xfine->Update(1.0,*xcoarse,0.0);
#if 0
     }
     else
     {
        cout << "**WRN** Maps are not equal in\n"
             << "**WRN** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        // create exporter from xcoarse->Map() to finegraph->RowMap()
        Epetra_Export* exporter = new Epetra_Export(xcoarse->Map(),xfine->Map());
        ierr = xfine->Export(*xcoarse,*exporter,Insert);
        if (ierr)
        {
           cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeF:\n"
                << "**ERR**: export from xcoarse to xfine returned err=" << ierr <<"\n"
                << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
        }
        delete exporter; exporter = 0;
     }
#endif
     // call fine level interface
     err = fineinterface_.computeF(*xfine,*Ffine,fillFlag);
     if (xfine) delete xfine; xfine = 0;
  
     // create importer from Ffine to F
     Epetra_Import* importer  = new Epetra_Import(F.Map(),Ffine->Map()); 

     // import FVecfine to FVec
     ierr = F.Import(*Ffine,*importer,Insert);
     if (ierr)
     {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeF:\n"
             << "**ERR**: import from Ffine to F returned err=" << ierr << "\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
     }

     //tidy up
     if (importer) delete importer; importer = 0;
     if (xcoarse)  delete xcoarse;  xcoarse  = 0;
     if (Ffine)    delete Ffine;    Ffine    = 0;
  }
  else // level_ > 0
  {
     // create Ffine and xfine matching the fine interface
     const Epetra_CrsGraph* finegraph = fineinterface_.getGraph();
     Epetra_Vector*         Ffine     = new Epetra_Vector(finegraph->RowMap(),false);
     Epetra_Vector*         xfine     = prolong_this_to_fine(*xcoarse);
     
     // call the fine grid user interface
     err = fineinterface_.computeF(*xfine,*Ffine,fillFlag);
     if (err==false)
     {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeF:\n"
             << "**ERR**: call to fine-userinterface returned false on level " << level_ << "\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
     }
     if (xfine) delete xfine; xfine = 0;
     
     // get the answer Ffine back to this coarse level
     Epetra_Vector* Fcoarse = restrict_fine_to_this(*Ffine);
     if (Ffine) delete Ffine; Ffine = 0;
     
     // FIXME:: after intensive testing, this test might not be necessary
#if 0
     samemap = F.Map().PointSameAs(Fcoarse->Map());
     if (samemap)
     {
#endif
        F.Update(1.0,*Fcoarse,0.0);
#if 0
     }
     else
     {
        cout << "**WRN** Maps are not equal in\n"
             << "**WRN** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        // create importer from Fcoarse to F 
        Epetra_Import* importer = new Epetra_Import(F.Map(),Fcoarse->Map()); 
        // import F from Fcoarse
        ierr = F.Import(*Fcoarse,*importer,Insert); 
        if (ierr)
        {
           cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeF:\n"
                << "**ERR**: import from Fcoarse to F returned err=" << ierr <<"\n"
                << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
        }
        if (importer) delete importer; importer = 0;
     }
#endif     
     if (Fcoarse) delete Fcoarse; Fcoarse = 0;
     if (xcoarse) delete xcoarse; xcoarse = 0;
  } // level_ > 0

  // check for FAS option
  if (isFAS()==true)
  {
     if (fbar_==NULL || fxbar_==NULL)
     {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeF:\n"
             << "**ERR**: isFAS is true and f-vector is NULL\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
     }
#if 0
     if (F.Map().SameAs(fbar_->Map())  != true ||
         F.Map().SameAs(fxbar_->Map()) != true)
     {
        cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeF:\n"
             << "**ERR**: mismatch in Maps of F and fbar_/fxbar_\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
     }
#endif
     F.Update(-1.0,*fxbar_,1.0,*fbar_,1.0);

  }

  return err;
}

/*----------------------------------------------------------------------*
 |  evaluate Jacobian           (public, derived)            m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Nox_CoarseProblem_Interface::computeJacobian(const Epetra_Vector& x)
{
  cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computeJacobian(...):\n"
       << "**ERR**: this  is NOT supposed to be called????????\n"
       << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  return true;
}

/*----------------------------------------------------------------------*
 |  compute user supplied preconditioner (public, derived)   m.gee 12/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::Nox_CoarseProblem_Interface::computePreconditioner(
                                                    const Epetra_Vector& x,
		                                    NOX::Parameter::List* precParams)
{
  cout << "**ERR**: ML_Epetra::Nox_CoarseProblem_Interface::computePreconditioner(...):\n"
       << "**ERR**: this  is NOT supposed to be called????????\n"
       << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  return(true);
}
//-----------------------------------------------------------------------------

#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
