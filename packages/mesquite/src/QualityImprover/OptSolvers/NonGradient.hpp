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
/*!
  \file   NonGradient.hpp
  \brief  

  The NonGradient Class implements the steepest descent algorithm in
  order to move a free vertex to an optimal position given an
  ObjectiveFunction object and a QualityMetric object.

*/

#ifndef Mesquite_NonGradient_hpp 
#define Mesquite_NonGradient_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "PatchSetUser.hpp"

namespace MESQUITE_NS
{
  class ObjectiveFunction;

  /*! \class NonGradient

      This is a basic implementation of the steepest descent optimization algorithm.
      It works on patches of any size.  The step size is hard-wired.
      This is only for testing purposes. */ 
  class NonGradient : public VertexMover, public PatchSetUser
  {
  public:
    MESQUITE_EXPORT 
    NonGradient(ObjectiveFunction* of);

    MESQUITE_EXPORT 
    NonGradient(ObjectiveFunction* of, MsqError& err );

    MESQUITE_EXPORT virtual
    ~NonGradient() { }
    
    MESQUITE_EXPORT virtual
    std::string get_name() const;
    
    MESQUITE_EXPORT virtual
    PatchSet* get_patch_set();
    
    MESQUITE_EXPORT
    bool project_gradient() const 
      { return projectGradient; }
    
    MESQUITE_EXPORT
    void project_gradient( bool yesno ) 
      { projectGradient = yesno; }

    int getDimension()
    { 
      return(mDimension);
    }
    double getThreshold()
    { 
      return(mThreshold);
    }
    double getTolerance()
    { 
      return(mTolerance);
    }
    int getMaxNumEval()
    { 
      return(mMaxNumEval);
    }
    void setDimension(int dimension)
    { 
      mDimension = dimension;
    }
    void setThreshold(double threshold)
    { 
      mThreshold = threshold;
    }
    void setTolerance(double ftol)
    { 
      mTolerance = ftol;
    }
    void setMaxNumEval(int maxNumEval)
    { 
      mMaxNumEval = maxNumEval;
    }
    void getRowSum( int numRow, int numCol, std::vector<double>& matrix, std::vector<double>& rowSum);
    bool testRowSum( int numRow, int numCol, double* matrix, double* rowSum);
    double evaluate( double localArray[], PatchData &pd, MsqError &err );
    //! edgeLenght is a length scale for the initial polytope.
    int initSimplex(double edgeLength); 
    //! matrix stored by column as a std::vector
    std::vector<double> simplex; 
    std::vector<double> height; 
    void printPatch( const PatchData &pd, MsqError &err );
    //! Obtain diagnostic data during optimization
    //! off=level 0, ... level 3 = maximal
    MESQUITE_EXPORT void set_debugging_level(int level)
    {
      mNonGradDebug=level;
    }

  protected:
    MESQUITE_EXPORT virtual
    void initialize( PatchData &pd, MsqError &err );
    MESQUITE_EXPORT virtual
    void optimize_vertex_positions( PatchData &pd, MsqError &err );
    MESQUITE_EXPORT virtual
    void initialize_mesh_iteration( PatchData &pd, MsqError &err );
    MESQUITE_EXPORT virtual
    void terminate_mesh_iteration( PatchData &pd, MsqError &err );
    MESQUITE_EXPORT virtual
    void cleanup();
  private:
    bool projectGradient;
    int mDimension;
    double mThreshold;// stop if       2(heightMax-heightMin)
    double mTolerance;//          ---------------------------------- < mTolerance
    int mMaxNumEval;  //          |heightMax|+|heightMin|+mThreshold
                      //      or numEval >= mMaxNumEval
    double amotry(std::vector<double>&, std::vector<double>& , double* , int , double, PatchData&, MsqError &err );
    int mNonGradDebug;

    NonGradient(const NonGradient &pd); //disable copying
    NonGradient& operator=(const NonGradient &pd);  //disable assignment
  };
  
}

#endif
