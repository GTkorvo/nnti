
//@HEADER
// ************************************************************************
//
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Epetra_Util.h"
#include "Trilinos_Util.h"
#include <string>

#include "Trilinos_Util_MatrixGallery.h"
#include "Trilinos_Util_ShellOptions.h"

// ================================================ ====== ==== ==== == =
Trilinos_Util_MatrixGallery::Trilinos_Util_MatrixGallery(const string name, 
							 const Epetra_Comm & comm ) :
  name_(name), comm_(&comm)
{
  ZeroOutData();
  // verbosity level
  if( comm_->MyPID()==0 ) verbose_ = true;
  else verbose_ = false;
  // fix error message
  ErrorMsg = "ERROR [Trilinos_Util_MatrixGallery]: ";
  OutputMsg = "Trilinos_Util_MatrixGallery: ";
  
}

// ================================================ ====== ==== ==== == =
Trilinos_Util_MatrixGallery::Trilinos_Util_MatrixGallery(const string name, 
							 const Epetra_Map & map ) :
  name_(name), comm_(&(map.Comm()))
{
  ZeroOutData();
  // verbosity level
  if( comm_->MyPID()==0 ) verbose_ = true;
  else verbose_ = false;
  // fix error message
  ErrorMsg = "ERROR [Trilinos_Util_MatrixGallery]: ";
  OutputMsg = "Trilinos_Util_MatrixGallery: ";
  
  map_ = new Epetra_Map(map);
  NumGlobalElements_ = map_->NumGlobalElements();
  NumMyElements_ = map_->NumMyElements();
  MyGlobalElements_ = map_->MyGlobalElements( );
    
}
  
// ================================================ ====== ==== ==== == =
Trilinos_Util_MatrixGallery::~Trilinos_Util_MatrixGallery(void) 
{

  // linear problem
  if( LinearProblem_ != NULL ) delete LinearProblem_;
  if( VbrLinearProblem_ != NULL ) delete VbrLinearProblem_;

  // VBR data
  if( VbrMatrix_ != NULL ) delete VbrMatrix_;  
  if( VbrExactSolution_ != NULL ) delete VbrExactSolution_;
  if( VbrStartingSolution_ != NULL ) delete VbrStartingSolution_;
  if( VbrRhs_ != NULL ) delete VbrRhs_;
  if( BlockMap_ != NULL ) delete BlockMap_;

  // Crs data
  if( matrix_ != NULL ) delete matrix_;
  if( ExactSolution_ != NULL ) delete ExactSolution_;
  if( StartingSolution_ != NULL ) delete StartingSolution_;
  if( rhs_ != NULL ) delete rhs_;
  if( map_ != NULL ) delete map_;

  // vectors
  if( VectorA_ != NULL ) delete VectorA_;
  if( VectorB_ != NULL ) delete VectorB_;
  if( VectorC_ != NULL ) delete VectorC_;
  if( VectorD_ != NULL ) delete VectorD_;
  if( VectorE_ != NULL ) delete VectorE_;
  if( VectorF_ != NULL ) delete VectorF_;
  if( VectorG_ != NULL ) delete VectorG_;

  // put to default values
  ZeroOutData();
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::Set(const string parameter, const int value)
{

  if( parameter == "problem_size" ) {
    if( value <= 0 ) {
      cerr << ErrorMsg << "problem size must be greater than 1\n";
      return -1;
    }
    if( map_ != NULL ) {
      cerr << ErrorMsg << "map object already set. Continuing with\n"
	   << ErrorMsg << "problemSize = " << NumGlobalElements_ << endl;
      return -2;
    }
    NumGlobalElements_ = value;
    return 0;

  } else if( parameter == "nx" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "nx must be greater than 0\n";
      return -1;
    }

    nx_ = value;
    return 0;
      
  } else if( parameter == "ny" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "ny must be greater than 0\n";
      return -1;
    }

    ny_ = value;
    return 0;
      
  } else if( parameter == "nz" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "nz must be greater than 0\n";
      return -1;
    }

    nz_ = value;
    return 0;
  } else if( parameter == "mx" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "mx must be greater than 0\n";
      return -1;
    }

    mx_ = value;
    return 0;
      
  } else if( parameter == "my" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "my must be greater than 0\n";
      return -1;
    }

    my_ = value;
    return 0;
      
  } else if( parameter == "mz" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "mz must be greater than 0\n";
      return -1;
    }

    mz_ = value;
    return 0;

  } else if( parameter == "num_pde_eqns" ) {

    if( value <= 0 ) {
      cerr << ErrorMsg << "num pde eqns must be greater than 0\n";
      return -1;
    }

    NumPDEEqns_ = value;
    return 0;
  } 

  cerr << ErrorMsg << "input string (" << parameter << ") not valid\n";
  return -2;

}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::Set(const string parameter, const string value )
{

  if( parameter == "problem_type" ) {
    name_ = value;
  }
  else if( parameter == "map_type" ) {
    MapType_ = value;
  }
  else if( parameter == "exact_solution" ) {
    ExactSolutionType_ = value;
  }
  else if( parameter == "matrix_name" ) {
    FileName_ = value;
  }    
  else if( parameter == "starting_solution" ) {
    StartingSolutionType_ = value;
  }
  else if( parameter == "output" ) {
    if( value == "none" ) verbose_ = false;
    if( value == "proc 0" ) {
      if( comm_->MyPID()==0 ) verbose_ = true;
      else verbose_ = false;
    } else {
      verbose_ = true;
    }
  } else if( parameter == "expand_type" ) {
    ExpandType_ = value;
  } else {
    cerr << ErrorMsg << "wrong input parameter (" << parameter << ")\n";
    return -1;
  }

  return 0;
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::Set(const string parameter, const double value)
{

  if( parameter == "a" ) {
    a_ = value;
    return 0;
  } else if( parameter == "b" ) {
    b_ = value;
    return 0;
  } else if( parameter == "c" ) {
    c_ = value;
    return 0;
  } else if( parameter == "d" ) {
    d_ = value;
    return 0;
  } else if( parameter == "e" ) {
    e_ = value;
    return 0;
  } else if( parameter == "f" ) {
    f_ = value;
    return 0;
  } else if( parameter == "g" ) {
    g_ = value;
    return 0;
  }

  cerr << ErrorMsg << "input string not valid\n";
  return -2;
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::Set(const string parameter, const Epetra_Vector & value)
{

  if( value.Map().SameAs(*map_) == false ) {
    cerr << ErrorMsg << "input vector must have the same map used to\n"
	 << ErrorMsg << "create the Trilinos_Util_MatrixGallery object. Continuing\n";
    return -2;
  }
    
  if( parameter == "a" ) {
    VectorA_ = new Epetra_Vector(value);
  } else if( parameter == "b" ) {
    VectorB_ = new Epetra_Vector(value);
  }
  else if( parameter == "c" ) {
    VectorC_ = new Epetra_Vector(value);
  }
  else if( parameter == "d" ) {
    VectorD_ = new Epetra_Vector(value);
  }
  else if( parameter == "e" ) {
    VectorE_ = new Epetra_Vector(value);
  }
  else if( parameter == "f" ) {
    VectorF_ = new Epetra_Vector(value);
  }
  else if( parameter == "g" ) {
    VectorG_ = new Epetra_Vector(value);
  } else {
    cerr << ErrorMsg << "input string not valid\n";
    return -3;
  }

  return 0;
}

// ================================================ ====== ==== ==== == =  

int Trilinos_Util_MatrixGallery::Set(Trilinos_Util_ShellOptions & S)
{
  // this can be done with STL: iterators and sets.
  // However, on some Sandia computers, STL is not that simple
  // to use. For this reason, I prefered this ugliest way

  int count;
  
  string Options[10];
  
  // all options with strings
  count = 0;
  Options[count++] = "problem_type";
  Options[count++] = "map_type";
  Options[count++] = "exact_solution";
  Options[count++] = "matrix_name";
  Options[count++] = "starting_solution";
  Options[count++] = "output";
  Options[count++] = "expand_type";
  
  for( int i=0 ; i<count ; i++ ) {
    string parameter = "-"+Options[i];    
    if( S.Have(parameter) == true ) {
      string value = S.GetString(parameter);
      Set(Options[i],value);
      
    }
  }

  // all options with integers
  Options[0] = "problem_size";
  Options[1] = "nx";
  Options[2] = "ny";
  Options[3] = "nz";
  Options[4] = "mx";
  Options[5] = "my";
  Options[6] = "mz";
  Options[7] = "num_pde_eqns";

  for(  int i=0 ; i<8 ; i++ ) {
    string parameter = "-"+Options[i];   
    if( S.Have(parameter) == true ) {
      Set(Options[i],S.GetInt(parameter));
    }
  }
  
  // all options with doubles
  Options[0] = "a";
  Options[1] = "b";
  Options[2] = "c";
  Options[3] = "d";
  Options[4] = "e";
  Options[5] = "f";
  Options[6] = "g";
  for( int i=0 ; i<7 ; i++ ) {
    string parameter = "-"+Options[i];   
    if( S.Have(parameter) == true ) {
      Set(Options[i],S.GetDouble(parameter));
    }
  }

  return 0;
}

// ================================================ ====== ==== ==== == =  
Epetra_CrsMatrix * Trilinos_Util_MatrixGallery::GetMatrix(void) 
{
  if( matrix_ == NULL ) assert(CreateMatrix()==0);
  return( matrix_ );
}

// ================================================ ====== ==== ==== == =  
Epetra_CrsMatrix & Trilinos_Util_MatrixGallery::GetMatrixRef(void) 
{
  if( matrix_ == NULL ) assert(CreateMatrix()==0);
  return( *matrix_ );
}
  
// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util_MatrixGallery::GetExactSolution(void) 
{
  if( ExactSolution_ == NULL ) assert(CreateExactSolution()==0);
  return ExactSolution_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util_MatrixGallery::GetStartingSolution(void)
{
  if( StartingSolution_ == NULL ) assert(CreateStartingSolution()==0);
  return StartingSolution_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util_MatrixGallery::GetRHS(void)
{
  if( rhs_ == NULL ) assert(CreateRHS()==0);
  return rhs_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util_MatrixGallery::GetVbrRHS(void)
{
  if( VbrRhs_ == NULL ) assert(CreateVbrRHS()==0);
  return VbrRhs_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util_MatrixGallery::GetVbrExactSolution(void)
{
  if( VbrExactSolution_ == NULL ) assert(CreateVbrExactSolution()==0);
  return VbrExactSolution_;
}

// ================================================ ====== ==== ==== == =
Epetra_Vector * Trilinos_Util_MatrixGallery::GetVbrStartingSolution(void)
{
  if( VbrStartingSolution_ == NULL ) assert(CreateVbrStartingSolution()==0);
  return VbrStartingSolution_;
}

// ================================================ ====== ==== ==== == =
const Epetra_Map * Trilinos_Util_MatrixGallery::GetMap(void)
{
  if( map_ == NULL ) assert(CreateMap()==0);
    
  return map_;
}

// ================================================ ====== ==== ==== == =
const Epetra_Map & Trilinos_Util_MatrixGallery::GetMapRef(void)
{
  if( map_ == NULL ) assert(CreateMap()==0);
    
  return *map_;
}

// ================================================ ====== ==== ==== == =
const Epetra_BlockMap * Trilinos_Util_MatrixGallery::GetBlockMap(void)
{
  if( BlockMap_ == NULL ) assert(CreateBlockMap()==0);
    
  return BlockMap_;
}

// ================================================ ====== ==== ==== == =
const Epetra_BlockMap & Trilinos_Util_MatrixGallery::GetBlockMapRef(void)
{
  if( BlockMap_ == NULL ) assert(CreateBlockMap()==0);
    
  return *BlockMap_;
}

// ================================================ ====== ==== ==== == =  
Epetra_VbrMatrix * Trilinos_Util_MatrixGallery::GetVbrMatrix(const int NumPDEEqns) 
{

  if( NumPDEEqns != NumPDEEqns_ ) {
    if( BlockMap_ != NULL ) {
      delete BlockMap_;
      BlockMap_ = NULL;
    }
    NumPDEEqns_ = NumPDEEqns;
      
  }

  return( GetVbrMatrix() );
    
}

// ================================================ ====== ==== ==== == =
Epetra_VbrMatrix * Trilinos_Util_MatrixGallery::GetVbrMatrix(void)
{
    
  if( VbrMatrix_ == NULL ) assert(CreateVbrMatrix()==0);

  return VbrMatrix_;
    
}

// ================================================ ====== ==== ==== == =
Epetra_VbrMatrix & Trilinos_Util_MatrixGallery::GetVbrMatrixRef(void)
{
    
  if( VbrMatrix_ == NULL ) assert(CreateVbrMatrix()==0);

  return *VbrMatrix_;
    
}

// ================================================ ====== ==== ==== == =
Epetra_LinearProblem * Trilinos_Util_MatrixGallery::GetLinearProblem(void) 
{
  // pointers, not really needed
  Epetra_CrsMatrix * A;
  Epetra_Vector * RHS;
  Epetra_Vector * StartingSolution;

  A = GetMatrix();
  RHS = GetRHS();
  StartingSolution = GetStartingSolution();

  // create linear problem
  if( LinearProblem_ != NULL ) delete LinearProblem_;
  LinearProblem_ = new Epetra_LinearProblem(A,StartingSolution,RHS);

  return LinearProblem_;

}

// ================================================ ====== ==== ==== == =
Epetra_LinearProblem * Trilinos_Util_MatrixGallery::GetVbrLinearProblem(void) 
{
  // pointers, not really needed
  Epetra_VbrMatrix * A;
  Epetra_Vector * RHS;
  Epetra_Vector * StartingSolution;

  A = GetVbrMatrix();
  RHS = GetVbrRHS();
  StartingSolution = GetVbrStartingSolution();

  // create linear problem
  if( VbrLinearProblem_ != NULL ) delete VbrLinearProblem_;
  VbrLinearProblem_ = new Epetra_LinearProblem(A,StartingSolution,RHS);

  return VbrLinearProblem_;

}
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::ComputeResidual(double & residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution_ are
  //  created by CreateRHS if necessary)
  if( rhs_ == NULL ) assert(CreateRHS()==0);

  Epetra_Vector Ax(*map_);
  assert(matrix_->Multiply(false, *StartingSolution_, Ax)==0);

  assert(Ax.Update(1.0, *rhs_, -1.0)==0);

  assert(Ax.Norm2(&residual)==0);

  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::ComputeDiffBetweenStartingAndExactSolutions(double & residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( rhs_ == NULL ) assert(CreateRHS()==0);

  Epetra_Vector temp(*map_);

  assert(temp.Update(1.0, *ExactSolution_, -1.0, *StartingSolution_, 0.0)==0);

  assert(temp.Norm2(&residual)==0);

  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMap(void)
{

  Epetra_Time Time(*comm_);

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating Map `" << MapType_ << "'...\n";
  }

  // first get the problem size. For some problems. the user can
  // specify the problem size using different parameters (e.g.,
  // nx and ny for a 2D Laplace problem). I need the internal
  // variable NumGlobalElements_ properly set before continuing.
  // NOTE: for HB problems, this value has already been set
  
  if( name_ == "diag" || name_ == "tridiag"  ||
      name_ == "laplace 1d" || name_ == "eye" ||
      name_ == "lehmer" || name_ == "minij" ||
      name_ == "ris" || name_ == "hilbert" ||
      name_ == "jordblock" ) {
    if( NumGlobalElements_ <= 0 ) {
      if( nx_ > 0 ) NumGlobalElements_ = nx_;
      else {
	cerr << ErrorMsg << "problem size not correct\n";
	return -1;
      }
    }
  }
    
  else if( name_ == "laplace_2d" || name_ == "cross_stencil_2d"
	   || name_ == "recirc_2d" ) {
    if( NumGlobalElements_ <= 0 ) {  
      if( nx_ > 0 && ny_ > 0 )
	NumGlobalElements_ = nx_*ny_;
      else {
	cerr << ErrorMsg << "problem size not correct\n";
	return -1;
      }
    }
  }
    
  else if( name_ == "laplace_3d" || name_ == "cross_stencil_3d" ) {
    if( NumGlobalElements_ <= 0 ) {
      if( nx_ > 0 && ny_ > 0 && nz_ > 0 )
	NumGlobalElements_ = nx_*ny_*nz_;
      else {
	cerr << ErrorMsg << "problem size not correct\n";
	return -1;
      }
    }

  } else if( name_ == "hb" ) {
    // The global number of elements has been set in ReadHBMatrix
    if( NumGlobalElements_ <= 0 ) {
      cerr << ErrorMsg << "problem size not correct\n";
      return -1;
    }
    
  } else {

    cerr << ErrorMsg << "matrix name is incorrect or not set ("
	 << name_ << ")\n";
    exit( EXIT_FAILURE );

  }

  // check out whether one is using only one proc or not.
  // If yes, creation of map is straightforward. Then return.
  
  if( comm_->NumProc() == 1 ) {

    map_ = new Epetra_Map(NumGlobalElements_,0,*comm_);

  } else {

    // Here below more than one processor.
  
    if( MapType_ == "linear" ) {
      
      map_ = new Epetra_Map (NumGlobalElements_,0,*comm_);
      
    } else if( MapType_ == "box" ) {

      if( mx_ == -1 || my_ == -1 ) {
	mx_ = (int)sqrt((double)(comm_->NumProc()));
	my_ = mx_;
	if( mx_ * my_ != comm_->NumProc() ) {
	  cerr << ErrorMsg << "number of processes must be square number\n"
	       << ErrorMsg << "otherwise set mx and my\n";
	  return -1;
	}
      }
      
      if( nx_ == -1 || ny_ == -1 ) {
	nx_ = (int)sqrt((double)NumGlobalElements_);
	ny_ = nx_;
	if( nx_ * ny_ != NumGlobalElements_ ) {
	  cerr << ErrorMsg << "number of elements must be square number\n"
	       << ErrorMsg << "otherwise set mx and my\n";
	  return -1;
	}
      }
      
      // how to divide the axis
      
      int modx = (nx_+(nx_%mx_))/mx_;
      int mody = (ny_+(ny_%my_))/my_;
      
      int MyPID = comm_->MyPID(), startx, starty, endx, endy;
      int xpid = MyPID/mx_;
      int ypid = MyPID%my_;
      
      startx = xpid*modx;
      if( (xpid+1)*modx < nx_ ) endx = (xpid+1)*modx;
      else endx = nx_;
      starty = ypid*mody;
      if( (ypid+1)*mody < ny_ ) endy = (ypid+1)*mody;
      else endy = ny_;
      
      int NumMyElements = (endx-startx)*(endy-starty);
      int * MyGlobalElements = new int[NumMyElements];
      int count = 0;
      
      for( int i=startx ; i<endx ; ++i ) {
	for( int j=starty ; j<endy ; ++j ) {
	  MyGlobalElements[count++] = i+j*nx_;
	}
      }
      
      map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,0,*comm_);
      
      // I delete this guy to that this case is not different from the
      // others, and I don't have to clean up this mess while
      // destroying the object.
      
      delete [] MyGlobalElements;
      
    } else if( MapType_ == "interlaced" ) {
      
      // this is the first funky map. Nodes are assigned so that
      // node 0 is given to proc 0, node 1 to proc 1, and
      // node i to proc i%NumProcs. Probably not the best, but it
      // results in decompositions with lots of boundary nodes.
      
      int NumProcs = comm_->NumProc();
      int MyPID = comm_->MyPID();
      
      int NumMyElements = NumGlobalElements_/NumProcs;
      if( MyPID < NumGlobalElements_%NumProcs ) NumMyElements++;
      
      int count = 0;
      int * MyGlobalElements = new int[NumMyElements];
      
      for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
	if( i%NumProcs == MyPID ) 
	  MyGlobalElements[count++] = i;
      }
      
      if( count != NumMyElements ) {
	cerr << ErrorMsg << "something went wrong in CreateMap\n";
	cerr << ErrorMsg << "count = " << count << ", NumMyElements = "
	     << NumMyElements << endl;
	exit( EXIT_FAILURE);
      }
      
      map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,0,*comm_);
      delete [] MyGlobalElements;
      
    } else if( MapType_ == "random" ) {
      
      // this is even funkier. Random decomposition of nodes into procs.
      // It should result in a more ordered decomposition than "interlaced"
      // This is the idea: I create the map on proc 0, then I broadcast
      // it to all procs. This is not very efficient, but saves some
      // MPI calls.
      
      int * part = new int[NumGlobalElements_];
      
      if( comm_->MyPID() == 0 ) {
	Epetra_Util Util;
	
	for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
	  unsigned int r = Util.RandomInt();	
	  part[i] = r%(comm_->NumProc());
	}
      }
      
      comm_->Broadcast(part,NumGlobalElements_,0);
      
      // count the elements assigned to this proc
      int NumMyElements = 0;
      for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
	if( part[i] == comm_->MyPID() ) NumMyElements++;
      }
      
      // get the loc2global list
      int * MyGlobalElements = new int[NumMyElements];
      int count = 0;
      for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
	if( part[i] == comm_->MyPID() ) MyGlobalElements[count++] = i;
      }
      
      map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,
			     0,*comm_);
      
      delete [] MyGlobalElements;
      delete [] part;
      
    } else {
      
      cerr << ErrorMsg << "MapType has an incorrect value (" << MapType_ << ")\n";
      return -1;
      
    }
  }

  // local number of rows
  NumMyElements_ = map_->NumMyElements();
  // get update list
  MyGlobalElements_ = map_->MyGlobalElements( );

  if( verbose_ == true ) {
    cout << OutputMsg << "Time to create Map: "
	 << Time.ElapsedTime() << " (s)\n";
  }

  return 0;
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrix(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating Matrix...\n";
  }
  
  // HB matrices are different, as their dimension has to be read before.
  // Here the idea is to read the matrix on proc 0, then build the
  // map, then redistribute it linearly

  if( name_ == "hb" ) {

    Epetra_Time Time(*comm_);
    ReadHBMatrix();
    if( verbose_ == true ) {
      cout << OutputMsg << "Time to create matrix: "
	   << Time.ElapsedTime() << " (s)\n";
    }
    
  } else {
      
    if( map_ == NULL ) assert(CreateMap()==0);

    Epetra_Time Time(*comm_);
    
    if( name_ == "diag" ) assert(CreateMatrixDiag()==0);

    else if( name_ == "eye" ) assert(CreateEye()==0);
      
    else if( name_ == "tridiag" ) assert(CreateMatrixTriDiag()==0);
      
    else if( name_ == "laplace_1d" ) assert(CreateMatrixLaplace1d()==0);
      
    else if( name_ == "laplace_2d" ) assert(CreateMatrixLaplace2d()==0);

    else if( name_ == "recirc_2d" ) assert(CreateMatrixRecirc2d()==0);
      
    else if( name_ == "laplace_3d" ) assert(CreateMatrixLaplace3d()==0);
      
    else if( name_ == "cross_stencil_2d" ) assert(CreateMatrixCrossStencil2d()==0);
      
    else if( name_ == "cross_stencil_3d" ) assert(CreateMatrixCrossStencil3d()==0);

    else if( name_ == "lehmer" ) assert(CreateMatrixLehmer()==0);

    else if( name_ == "minij" ) assert(CreateMatrixMinij()==0);

    else if( name_ == "ris" ) assert(CreateMatrixRis()==0);

    else if( name_ == "hilbert" ) assert(CreateMatrixHilbert()==0);

    else if( name_ == "jordblock" ) assert(CreateMatrixJordblock()==0);
    
    else {
      cerr << ErrorMsg << "matrix name is incorrect or not set ("
	   << name_ << ")\n";
      exit( EXIT_FAILURE );
    }

    if( verbose_ == true ) {
      cout << OutputMsg << "Time to create matrix: "
	   << Time.ElapsedTime() << " (s)\n";
    }
  }

  return 0;    
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateVbrExactSolution(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating exact solution (VBR)...\n";
  }

  // check if already have one; in this case delete and rebuild
  if( VbrExactSolution_ != NULL ) delete VbrExactSolution_;
  // need an exact solution for Crs first
  if( ExactSolution_ == NULL ) assert(CreateExactSolution()==0);
  // need a block map first
  if( BlockMap_ == NULL ) assert(CreateBlockMap()==0);
  // now we can expand to the Vbr format
  VbrExactSolution_ = new Epetra_Vector(*BlockMap_);
  for( int j=0 ; j<NumMyElements_ ; j++ ) 
    for( int i=0 ; i<NumPDEEqns_ ; ++i ) {
      (*VbrExactSolution_)[j*NumPDEEqns_+i] = (*ExactSolution_)[j];
    }

  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateExactSolution(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating exact solution `"
	 << ExactSolutionType_ << "'...\n";
  }
  
  if( map_ == NULL ) assert(CreateMap()==0);

  if( ExactSolution_ == NULL ) {
    ExactSolution_ = new Epetra_Vector(*map_);
    if( ExactSolutionType_ == "random" ) {
      ExactSolution_->Random();
    } else if( ExactSolutionType_ == "constant" ) {
      ExactSolution_->PutScalar(1.0);
    } else if( ExactSolutionType_ == "linear" ) {
      for( int i=0 ; i<NumMyElements_ ; i++ ) 
	(*ExactSolution_)[i] = alpha_*MyGlobalElements_[i];
    } else {
      cerr << ErrorMsg << "exact solution type is not correct : "
	   << ExactSolutionType_ << endl;
      return -1;
    }
  }

  return 0;
  
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateStartingSolution(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating starting solution `"
	 << StartingSolutionType_ << "'...\n";
  }

  if( map_ == NULL ) assert(CreateMap()==0);

  if( StartingSolution_ == NULL ) {
    StartingSolution_ = new Epetra_Vector(*map_);
    if( StartingSolutionType_ == "random" ) {
      StartingSolution_->Random();
    } else if( StartingSolutionType_ == "zero" ) {
      StartingSolution_->PutScalar(0.0);
    } else {
      cerr << ErrorMsg << "starting solution type is not correct : "
	   << StartingSolutionType_ << endl;
      return -1;
    }
  }

  return 0;
  
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateVbrStartingSolution(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating Starting Solution (VBR)...\n";
  }

  if( VbrStartingSolution_ != NULL ) {
    delete VbrStartingSolution_;
    VbrStartingSolution_ = NULL;
  }
    
  // need a rhs for crs
  if( StartingSolution_ == NULL ) assert(CreateStartingSolution()==0);
  // need a block map based on map_
  if( BlockMap_ == NULL ) assert(CreateBlockMap()==0);
  // now we can expand to the Vbr format
  VbrStartingSolution_ = new Epetra_Vector(*BlockMap_);
  for( int j=0 ; j<NumMyElements_ ; j++ ) 
    for( int i=0 ; i<NumPDEEqns_ ; ++i ) {
      (*VbrStartingSolution_)[j*NumPDEEqns_+i] = (*StartingSolution_)[j];
    }

  return 0;
  
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateRHS(void)
{

  if( map_ == NULL ) assert(CreateMap()==0);
  if( matrix_ == NULL ) assert(CreateMatrix()==0);
  if( ExactSolution_ == NULL )  assert(CreateExactSolution()==0);

  if( rhs_ != NULL ) delete rhs_;
  
  Epetra_Time Time(*comm_);
    
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating RHS...\n";
  }
  
  rhs_ = new Epetra_Vector(*map_);
  matrix_->Multiply(false,*ExactSolution_,*rhs_);
    
  if( verbose_ == true ) {
    cout << OutputMsg << "Time to create RHS (matvec): "
         << Time.ElapsedTime() << " (s)\n";
  }

  return 0;
  
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateVbrRHS(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating RHS (VBR)...\n";
  }

  if( VbrRhs_ != NULL ) {
    delete VbrRhs_;
    VbrRhs_ = NULL;
  }
    
  // need a rhs for crs
  if( rhs_ == NULL ) assert(CreateRHS()==0);
  // need a block map based on map_
  if( BlockMap_ == NULL ) assert(CreateBlockMap()==0);
  // require VbrMatrix to be formed first
  if( VbrMatrix_ == NULL ) assert(CreateVbrMatrix()==0);
  // also need an exact solution
  if( VbrExactSolution_ == NULL )  assert(CreateVbrExactSolution()==0);

  VbrRhs_ = new Epetra_Vector( *BlockMap_);
  VbrMatrix_->Multiply(false,*VbrExactSolution_,*VbrRhs_);

  return 0;
  
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateEye(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `eye'...\n";
  }
  
  a_ = 1.0;
  assert(CreateMatrixDiag()==0);
  return 0;    
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixDiag(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `diag'...\n";
  }
  
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,1);
  double Value;
    
  for( int i=0 ; i<NumMyElements_; ++i ) {

    int Indices = MyGlobalElements_[i];
    if( VectorA_ == NULL ) Value = a_;
    else Value = (*VectorA_)[i];
      
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, &Value, &Indices)==0);
      
  }
    
  assert(matrix_->FillComplete()==0);

  return 0;
    
}
  
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixTriDiag(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `tridiag'...n";
  }
  
  int ierr;

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,3);

  double *Values = new double[2];
  Values[0] = b_; Values[1] = c_;
  int *Indices = new int[2];
  int NumEntries;

  for( int i=0 ; i<NumMyElements_; ++i ) {
    if (MyGlobalElements_[i]==0) {
      Indices[0] = 1;
      NumEntries = 1;
      if( VectorC_ == NULL ) Values[0] = c_;
      else Values[0] = (*VectorC_)[i];
    } else if (MyGlobalElements_[i] == NumGlobalElements_-1) {
      Indices[0] = NumGlobalElements_-2;
      NumEntries = 1;
      if( VectorC_ == NULL ) Values[0] = c_;
      else Values[0] = (*VectorC_)[i];
    } else {
      Indices[0] = MyGlobalElements_[i]-1;
      Indices[1] = MyGlobalElements_[i]+1;
      NumEntries = 2;
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, Values, Indices)==0);
    // Put in the diagonal entry
    if( VectorA_ == NULL ) Values[0] = a_;
    else Values[0] = (*VectorA_)[i];
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, Values, MyGlobalElements_+i)==0);
  }
  
  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  assert(matrix_->FillComplete()==0);

  delete [] Values;
  delete [] Indices;

  return 0;
    
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixLaplace1d(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `laplace_1d'...\n";
  }

  a_ = 2.0;
  b_ = -1.0;
  c_ = -1.0;

  assert(CreateMatrixTriDiag()==0);
  
  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixCrossStencil2d(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `cross_stencil_2d'...\n";
  }

  if( nx_ == -1 || ny_ == -1 ) {
    nx_ = (int)sqrt((double)NumGlobalElements_);
    ny_ = nx_;
      
    if( nx_ * ny_ != NumGlobalElements_ ) {
      cerr << ErrorMsg << "The number of global elements must be a square number\n"
	   << ErrorMsg << "(now is " << NumGlobalElements_ << "). Returning...\n";
      return -2;
    }
  }
    
  int left, right, lower, upper;
    
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,5);
    
  // Add  rows one-at-a-time
    
  double Values[4], diag;
  int Indices[4];
  int NumEntries;

  //    e
  //  b a c
  //    d
  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian2d(  MyGlobalElements_[i], nx_, ny_, 
			       left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      if( VectorB_ == NULL ) Values[NumEntries] = b_;
      else Values[NumEntries] = (*VectorB_)[i];
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      if( VectorC_ == NULL ) Values[NumEntries] = c_;
      else Values[NumEntries] = (*VectorC_)[i];
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      if( VectorD_ == NULL ) Values[NumEntries] = d_;
      else Values[NumEntries] = (*VectorD_)[i];
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      if( VectorE_ == NULL ) Values[NumEntries] = e_;
      else Values[NumEntries] = (*VectorE_)[i];
      ++NumEntries;
    }
    // put the off-diagonal entries
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, 
				       Values, Indices)==0);
    // Put in the diagonal entry
    if( VectorA_ == NULL ) diag = a_;
    else diag = (*VectorA_)[i];
	
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, 
				       &diag, MyGlobalElements_+i)==0);
  }
  matrix_->FillComplete();

  return 0;
      
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixLaplace2d(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `laplace_2d'...\n";
  }

  a_ = 4.0;
  b_ = -1.0;
  c_ = -1.0;
  d_ = -1.0;
  e_ = -1.0;

  assert(CreateMatrixCrossStencil2d()==0);

  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixRecirc2d(void)
{

  // largely based on a ML example by Ray Tuminaro
  
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `recirc_2d'...\n";
  }

  if( nx_ == -1 || ny_ == -1 ) {
    nx_ = (int)sqrt((double)NumGlobalElements_);
    ny_ = nx_;
      
    if( nx_ * ny_ != NumGlobalElements_ ) {
      cerr << ErrorMsg << "The number of global elements must be a square number\n"
	   << ErrorMsg << "(now is " << NumGlobalElements_ << "). Returning...\n";
      return -2;
    }
  }
  
  // need vectors for this job
  if( VectorA_ == NULL )  VectorA_ = new Epetra_Vector(*map_);
  if( VectorB_ == NULL )  VectorB_ = new Epetra_Vector(*map_);
  if( VectorC_ == NULL )  VectorC_ = new Epetra_Vector(*map_);
  if( VectorD_ == NULL )  VectorD_ = new Epetra_Vector(*map_);
  if( VectorE_ == NULL )  VectorE_ = new Epetra_Vector(*map_);

  VectorA_->PutScalar(0.0);
  VectorB_->PutScalar(0.0);
  VectorC_->PutScalar(0.0);
  VectorD_->PutScalar(0.0);
  VectorE_->PutScalar(0.0);
  
  double hx = 1.0/(nx_+1);
  double hy = 1.0/(ny_+1);

  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int ix, iy;
    ix = (MyGlobalElements_[i])%nx_;
    iy = (MyGlobalElements_[i] - ix)/nx_;
    double x = hx*(ix+1);
    double y = hy*(iy+1);
    double ConvX = a_*4*x*(x-1.)*(1.-2*y)/hx;
    double ConvY = -a_*4*y*(y-1.)*(1.-2*x)/hy;

    // convection part
    
    if( ConvX<0 ) {
      (*VectorE_)[i] += ConvX;
      (*VectorA_)[i] -= ConvX;
    } else {
      (*VectorD_)[i] -= ConvX;
      (*VectorA_)[i] += ConvX;
    }

    if( ConvY<0 ) {
      (*VectorB_)[i] += ConvY;
      (*VectorA_)[i] -= ConvY;
    } else {
      (*VectorC_)[i] -= ConvY;
      (*VectorA_)[i] += ConvY;
    }

    // add diffusion part
    (*VectorA_)[i] += b_*2./(hx*hx) + b_*2./(hy*hy);
    (*VectorB_)[i] -= b_/(hx*hx);
    (*VectorC_)[i] -= b_/(hx*hx);
    (*VectorD_)[i] -= b_/(hy*hy);
    (*VectorE_)[i] -= b_/(hy*hy);
      
    
  }
  
  assert(CreateMatrixCrossStencil2d()==0);

  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixLaplace3d(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `laplace_3d'...\n";
  }
  
  a_ = 6.0;
  b_ = -1.0;
  c_ = -1.0;
  d_ = -1.0;
  e_ = -1.0;
  f_ = -1.0;
  g_ = -1.0;

  assert(CreateMatrixCrossStencil3d()==0);

  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixCrossStencil3d(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `cross_stencil_3d'...\n";
  }
  
  if( nx_ == -1 || ny_ == -1 || nz_ == -1 ) {
    nx_ = (int)pow(1.0*NumGlobalElements_,0.333334);
    ny_ = nx_;
    nz_ = nx_;
      
    if( nx_ * ny_ *nz_ != NumGlobalElements_ ) {
      cerr << ErrorMsg << "The number of global elements must be a perfect cube\n"
	   << ErrorMsg << "(now is " << NumGlobalElements_ << "). Returning...\n";
      return -2;
    }
  }
    
  int left, right, lower, upper, below, above;
    
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,7);
    
  // Add  rows one-at-a-time
    
  double Values[6], diag;
  int Indices[6];
  int NumEntries;

  //    e 
  //  b a c
  //    d
  // + f below and g above
    
  for( int i=0 ; i<NumMyElements_; ++i ) {
    int NumEntries=0;
    GetNeighboursCartesian3d(MyGlobalElements_[i], nx_, ny_, nz_,
			     left, right, lower, upper, below, above);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      if( VectorB_ == NULL ) Values[NumEntries] = b_;
      else Values[NumEntries] = (*VectorB_)[i];
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      if( VectorC_ == NULL ) Values[NumEntries] = c_;
      else Values[NumEntries] = (*VectorC_)[i];
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      if( VectorD_ == NULL ) Values[NumEntries] = d_;
      else Values[NumEntries] = (*VectorD_)[i];
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      if( VectorE_ == NULL ) Values[NumEntries] = e_;
      else Values[NumEntries] = (*VectorE_)[i];
      ++NumEntries;
    }
    if( below != -1 ) {
      Indices[NumEntries] = below;
      if( VectorF_ == NULL ) Values[NumEntries] = f_;
      else Values[NumEntries] = (*VectorF_)[i];
      ++NumEntries;
    }
    if( above != -1 ) {
      Indices[NumEntries] = above;
      if( VectorG_ == NULL ) Values[NumEntries] = g_;
      else Values[NumEntries] = (*VectorG_)[i];
      ++NumEntries;
    }
    // put the off-diagonal entries
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], NumEntries, 
				       Values, Indices)==0);
    // Put in the diagonal entry
    if( VectorA_ == NULL ) diag = a_;
    else diag = (*VectorA_)[i];
	
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i], 1, 
				       &diag, MyGlobalElements_+i)==0);
  }

  matrix_->FillComplete();
  return 0;
      
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixLehmer(void)
{

    if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `lehmer'...\n";
  }

  // this is actually a dense matrix, stored into Crs format
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  for( int i=0 ; i<NumGlobalElements_ ; ++i ) Indices[i] = i;
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      if( i>=j ) Values[j] = 1.0*(j+1)/(i+1);
      else       Values[j] = 1.0*(i+1)/(j+1);
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumGlobalElements_, Values, Indices)==0);
      
  }
    
  delete [] Indices;
  delete [] Values;

  assert(matrix_->FillComplete()==0);
  
  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixMinij(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `minij'...\n";
  }
  
  // this is actually a dense matrix, stored into Crs format
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  for( int i=0 ; i<NumGlobalElements_ ; ++i ) Indices[i] = i;
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      if( i>=j ) Values[j] = 1.0*(j+1);
      else       Values[j] = 1.0*(i+1);
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumGlobalElements_, Values, Indices)==0);
      
  }
    
  delete [] Indices;
  delete [] Values;

  assert(matrix_->FillComplete()==0);
  
  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixRis(void)
{
  
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `ris'...\n";
  }

  // this is actually a dense matrix, stored into Crs format
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  for( int i=0 ; i<NumGlobalElements_ ; ++i ) Indices[i] = i;
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      Values[j] = 0.5/(NumGlobalElements_ -(i+1)-(j+1)+1.5);
      
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumGlobalElements_, Values, Indices)==0);
      
  }
    
  delete [] Indices;
  delete [] Values;

  assert(matrix_->FillComplete()==0);
  
  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixHilbert(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `hilbert'...\n";
  }

  // this is actually a dense matrix, stored into Crs format
  matrix_ = new Epetra_CrsMatrix(Copy,*map_,NumGlobalElements_);

  int * Indices = new int[NumGlobalElements_];
  double * Values = new double[NumGlobalElements_];

  for( int i=0 ; i<NumGlobalElements_ ; ++i ) Indices[i] = i;
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      Values[j] = 1.0/((i+1)+(j+1)-1);
    }
    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumGlobalElements_, Values, Indices)==0);
      
  }
    
  delete [] Indices;
  delete [] Values;

  assert(matrix_->FillComplete()==0);
  
  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateMatrixJordblock(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating matrix `jordblock'...\n";
  }

  matrix_ = new Epetra_CrsMatrix(Copy,*map_,2);

  int Indices[2];
  double Values[2];
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    int NumEntries = 0;
    if( MyGlobalElements_[i] != NumGlobalElements_-1 ) {
      Indices[NumEntries] = MyGlobalElements_[i]+1;
      Values[NumEntries] = 1.0;
      NumEntries++;
    }
    // diagonal contribution
    Indices[NumEntries] = MyGlobalElements_[i];
    if( VectorA_ != NULL ) Values[NumEntries] = (*VectorA_)[i];
    else                   Values[NumEntries] = a_;
    NumEntries++;

    assert(matrix_->InsertGlobalValues(MyGlobalElements_[i],
				       NumEntries, Values, Indices)==0);
      
  }
    
  assert(matrix_->FillComplete()==0);
  
  return 0;
}
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::ReadHBMatrix(void)
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Reading HB matrix `"
	 << FileName_ << "'...\n";
  }

  Epetra_Map * readMap;
  Epetra_CrsMatrix * readA; 
  Epetra_Vector * readx; 
  Epetra_Vector * readb;
  Epetra_Vector * readxexact;
    
  // Call routine to read in HB problem
    
  Trilinos_Util_ReadHb2Epetra((char*)FileName_.c_str(), *comm_, readMap, readA, readx, 
			      readb, readxexact);
    
  NumGlobalElements_ = readMap->NumGlobalElements();

  if( map_ != NULL ) delete map_;
  
  // create map for matrix. Use the normal function CreateMap
  // if the user has not specified "greedy" as map type.
  // In this latter case, form on proc 0 a map corresponding to
  // the greedy algorithm. This is some kind of graph decomposition
  // stuff, only cheaper (and that does not required any external
  // library)

  if( MapType_ == "greedy" ) {

    int * part = new int[NumGlobalElements_];

    if( comm_->MyPID() == 0 ) {

      int NumProcs = comm_->NumProc();
      int * ElementsPerDomain = new int[NumProcs];
      int * count = new int[NumProcs];

      // define how many nodes have to be put on each proc
      
      int div = NumGlobalElements_/NumProcs;
      int mod = NumGlobalElements_%NumProcs;

      for( int i=0 ; i<NumProcs ; ++i ) {
	count[i] = 0;
	ElementsPerDomain[i] = div;
	if( i<mod ) ElementsPerDomain[i]++;
      }
      
      for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
	part[i] = -1;
      }
      
      int MaxNnzPerRow = readA->MaxNumEntries();
      if( MaxNnzPerRow == 0 ) {
	cerr << ErrorMsg << "something went wrong in `CreateMatrix'\n"
	     << ErrorMsg << "MaxNnzPerRow == 0 \n";
	exit( EXIT_FAILURE );
      }

      int CrsNumEntries;
      int * CrsIndices;
      double * CrsValues;

      // start from row 0, assigned to domain 0
      int RootNode = 0;
      part[0] = 0;      
      int CurrentDomain = 0;
      
      bool ok = true;
      
      while( ok == true ) {

	readA->ExtractMyRowView(RootNode,CrsNumEntries,
				CrsValues,CrsIndices);

	ok = false;
	
	for( int j=0 ; j<CrsNumEntries ; ++j ) {

	  if( count[CurrentDomain] == ElementsPerDomain[CurrentDomain] ) {
	    CurrentDomain++;
	  }
	  
	  if( part[CrsIndices[j]] == -1 ) {
	    part[CrsIndices[j]] = CurrentDomain;
	    if( ok == false ) {
	      ok = true;
	      RootNode = CrsIndices[j];
	    }
	    count[CurrentDomain]++;
	  }
	}

	// check if some -1 nodes are still available
	if( ok == false ) {
	  for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
	    if( part[j] == -1 ) {
	      RootNode = j;
	      ok = true;
	      break;
	    }
	  }
	}
	      
      }

      delete [] ElementsPerDomain;
      delete [] count;
	
    }

    // now broadcast on all procs. This might be pretty expensive...
    comm_->Broadcast(part,NumGlobalElements_,0);

    for( int j=0 ; j<NumGlobalElements_ ; ++j ) {
      if( part[j] == -1 ) {
	cerr << ErrorMsg << "part[" << j << "] = -1 \n";
      }
    }    

    // count the elements assigned to this proc
    int NumMyElements = 0;
    for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
      if( part[i] == comm_->MyPID() ) NumMyElements++;
    }

    // get the loc2global list
    int * MyGlobalElements = new int[NumMyElements];
    int count = 0;
    for( int i=0 ; i<NumGlobalElements_ ; ++i ) {
      if( part[i] == comm_->MyPID() ) MyGlobalElements[count++] = i;
    }

    map_ = new Epetra_Map (NumGlobalElements_,NumMyElements,MyGlobalElements,
			   0,*comm_);
    
    delete [] MyGlobalElements;
    delete [] part;

  } else {    
    assert(CreateMap()==0);
  }

  // Create Exporter to distribute read-in matrix and vectors
  Epetra_Export exporter(*readMap, *map_);
  matrix_ = new Epetra_CrsMatrix(Copy, *map_, 0);
  StartingSolution_ = new Epetra_Vector(*map_);
  rhs_ = new Epetra_Vector(*map_);
  ExactSolution_ = new Epetra_Vector(*map_);
  StartingSolution_->Export(*readx, exporter, Add);
  rhs_->Export(*readb, exporter, Add);
  ExactSolution_->Export(*readxexact, exporter, Add);
  matrix_->Export(*readA, exporter, Add);
  
  assert(matrix_->FillComplete()==0);    

  delete readA;
  delete readx;
  delete readb;
  delete readxexact;
  delete readMap;
  
  // local number of rows
  NumMyElements_ = map_->NumMyElements();
  // get update list
  MyGlobalElements_ = map_->MyGlobalElements( );
  
  return 0;

}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateBlockMap(void) 
{
        
  if( verbose_ == true ) {
    cout << OutputMsg << "Creating BlockMap...\n";
  }

  if( map_ == NULL ) assert(CreateMap()==0);

  Epetra_Time Time(*comm_);
    
  if( NumPDEEqns_ <= 0 ) {
    cerr << ErrorMsg << "NumPDEEqns not correct (" << NumPDEEqns_ << "(\n";
    cerr << ErrorMsg << "Set it to 1\n";
    NumPDEEqns_ = 1;
  }

  MaxBlkSize_ = NumPDEEqns_;
  
  BlockMap_ = new Epetra_BlockMap(NumGlobalElements_,NumMyElements_,
				  MyGlobalElements_, 
				  NumPDEEqns_,0,*comm_);

  if( verbose_ == true ) {
    cout << OutputMsg << "Time to create BlockMap: "
	 << Time.ElapsedTime() << " (s)\n";
  }

  return 0;
  
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::CreateVbrMatrix(void) 
{

  if( verbose_ == true ) {
    cout << OutputMsg << "Creating VBR matrix...\n";
  }

  if( matrix_ == NULL ) assert(CreateMatrix()==0);
  if( BlockMap_ == NULL ) assert(CreateBlockMap()==0);

  int MaxNnzPerRow = matrix_->MaxNumEntries();
  if( MaxNnzPerRow == 0 ) {
    cerr << ErrorMsg << "something went wrong in `CreateMatrix'\n"
	 << ErrorMsg << "MaxNnzPerRow == 0 \n";
    return -3;
  }
    
  // create a VBR matrix based on BlockMap
  VbrMatrix_ = new Epetra_VbrMatrix(Copy, *BlockMap_,MaxNnzPerRow);

  // size of each VBR block
  int MaxBlockSize = MaxBlkSize_*MaxBlkSize_;

  int CrsNumEntries;
  int * CrsIndices;
  double * CrsValues;
    
  int * VbrIndices = new int[MaxNnzPerRow];
  double * VbrValues = new double[MaxBlockSize];
  int BlockRows = NumPDEEqns_;
  int ierr;
    
  // cycle over all the local rows. 
  
  for( int i=0 ; i<NumMyElements_ ; ++i ) {
    
    // get GID of local row
    int GlobalNode = MyGlobalElements_[i];
    // extract Crs row

    ierr = matrix_->ExtractMyRowView(i,CrsNumEntries,
				     CrsValues,CrsIndices);

    // matrix_ is in local form. Need global indices
    for( int kk=0 ; kk<CrsNumEntries ; ++kk) 
      VbrIndices[kk] = matrix_->GCID(CrsIndices[kk]);
    
    // with VBR matrices, we have to insert one block at time.
    // This required two more instructions, one to start this
    // process (BeginInsertGlobalValues), and another one to
    // commit the end of submissions (EndSubmitEntries).
    
    VbrMatrix_->BeginInsertGlobalValues(GlobalNode, CrsNumEntries, VbrIndices);

    int ExpandTypeInt;
    
    if( ExpandType_ == "zero_off_diagonal" ) ExpandTypeInt=0;
    else if( ExpandType_ == "random_off_diagonal" ) ExpandTypeInt=1;
    else {
      cerr << ErrorMsg << "ExpandType not correct (" << ExpandType_ << "\n";
      exit( EXIT_FAILURE );
    }
    Epetra_Util Util;
    double r;
    
    for( int i=0 ; i<CrsNumEntries ; ++i ) {
	
      for( int k=0 ; k<BlockRows ; ++k ) { // rows
	for( int h=0 ; h<BlockRows ; ++h ) { // cols
	  if( k == h ) VbrValues[k+h*BlockRows] = CrsValues[i];
	  else {
	    switch( ExpandTypeInt ) {
	    case 0:
	      r = 0.0;
	      break;
	    case 1:
	      // get a double between -1 and 1
	      r = Util.RandomDouble();
	      // scale it so that the sum of the block off-diagonal
	      // is not greater than the block diangonal
	      r /= (1.5*CrsValues[i]*BlockRows);
	      break;
	    }
	    VbrValues[k+h*BlockRows] = r; 
	  }
	}
      }
	  
      VbrMatrix_->SubmitBlockEntry(VbrValues,BlockRows,BlockRows,BlockRows);
	
    }

    VbrMatrix_->EndSubmitEntries();
  }
    
  delete [] VbrIndices;
  delete [] VbrValues;

  VbrMatrix_->FillComplete();

  return 0;
}

// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::ComputeResidualVbr(double & residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( VbrRhs_ == NULL ) assert(CreateVbrRHS()==0);

  Epetra_Vector Ax(*BlockMap_);
  assert(VbrMatrix_->Multiply(false, *VbrStartingSolution_, Ax)==0);

  assert(Ax.Update(1.0, *VbrRhs_, -1.0)==0);

  assert(Ax.Norm2(&residual)==0);

  return 0;
}
// ================================================ ====== ==== ==== == =
int Trilinos_Util_MatrixGallery::ComputeDiffBetweenStartingAndExactSolutionsVbr(double & residual)
{

  // create solution and rhs if needed (matrix_ and ExactSolution are
  //  created by CreateRHS is necessary)
  if( VbrRhs_ == NULL ) assert(CreateVbrRHS()==0);

  Epetra_Vector temp(*BlockMap_);

  assert(temp.Update(1.0, *VbrExactSolution_, -1.0, *VbrStartingSolution_, 0.0)==0);

  assert(temp.Norm2(&residual)==0);

  return 0;
}

// ================================================ ====== ==== ==== == =  
void Trilinos_Util_MatrixGallery::GetNeighboursCartesian2d( const int i, const int nx, const int ny,
				int & left, int & right, 
				int & lower, int & upper) 
{

  int ix, iy;
  ix = i%nx;
  iy = (i - ix)/nx;

  if( ix == 0 ) 
    left = -1;
  else 
    left = i-1;
  if( ix == nx-1 ) 
    right = -1;
  else
    right = i+1;
  if( iy == 0 ) 
    lower = -1;
  else
    lower = i-nx;
  if( iy == ny-1 ) 
    upper = -1;
  else
    upper = i+nx;

  return;

}

// ================================================ ====== ==== ==== == =
void Trilinos_Util_MatrixGallery::GetNeighboursCartesian3d( const int i, const int nx, const int ny, const int nz,
				int & left, int & right, int & lower, int & upper,
				int & below, int & above ) 
{

  int ixy, iz;
  ixy = i%(nx*ny);
    
  iz = (i - ixy)/(nx*ny);

  if( iz == 0 ) 
    below = -1;
  else 
    below = i-nx*ny;
  if( iz == nz-1 ) 
    above = -1;
  else
    above = i+nx*ny;

  GetNeighboursCartesian2d( ixy, nx, ny, left, right, lower, upper);
    
  if( left != -1 ) left += iz*(nx*ny);
  if( right != -1 ) right += iz*(nx*ny);
  if( lower != -1 ) lower += iz*(nx*ny);
  if( upper != -1 ) upper += iz*(nx*ny);

  return;

}

// ================================================ ====== ==== ==== == =
void Trilinos_Util_MatrixGallery::ZeroOutData() 
{
  NumGlobalElements_ = -1;
  nx_ = -1;    ny_ = -1;     nz_ = -1;
  mx_ = -1;    mx_ = -1;     mz_ = -1;
    
  a_ = 1.0, b_ = 0.0, c_ = 0.0, d_ = 0.0, e_ = 0.0, f_ = 0.0, g_ = 0.0;
  alpha_ = 1.0;
  beta_ = 0.0;
  gamma_ = 0.0;
  delta_ = 0.0;
    
  VectorA_ = NULL;
  VectorB_ = NULL;
  VectorC_ = NULL;
  VectorD_ = NULL;
  VectorE_ = NULL;
  VectorF_ = NULL;
  VectorG_ = NULL;
    
  map_ = NULL;
  matrix_ = NULL;
  ExactSolution_ = NULL;
  StartingSolution_ = NULL;
  rhs_ = NULL;

  BlockMap_ = NULL;
  VbrMatrix_ = NULL;
  VbrExactSolution_ = NULL;
  VbrStartingSolution_ = NULL;
  VbrRhs_ = NULL;
  
  MapType_ = "linear";
  ExactSolutionType_ = "constant";
  StartingSolutionType_ = "zero";
  ExpandType_ = "zero_off_diagonal";

  NumPDEEqns_= 1;

  LinearProblem_ = NULL;
  VbrLinearProblem_ = NULL;
    
}

void Trilinos_Util_MatrixGallery::PrintMatrixAndVectors() 
{
  PrintMatrixAndVectors(cout);
}

void Trilinos_Util_MatrixGallery::PrintMatrixAndVectors(ostream & os) 
{

  if( comm_->MyPID() == 0 ) {
    os << "*** MATRIX ***\n";
  }

  os << *matrix_;

  if( comm_->MyPID() == 0 ) {
    os << "*** RHS ***\n";
  }
  
  os << *rhs_;

  return;

}

void Trilinos_Util_MatrixGallery::PrintVbrMatrixAndVectors(ostream & os)
{

  if( comm_->MyPID() == 0 ) {
    os << "*** MATRIX (VBR) ***\n";
  }

  os << *VbrMatrix_;

  if( comm_->MyPID() == 0 ) {
    os << "*** RHS (VBR) ***\n";
  }

  os << *VbrRhs_;

  return;

}

void Trilinos_Util_MatrixGallery::PrintVbrMatrixAndVectors() 
{
  PrintVbrMatrixAndVectors(cout);
}

ostream & operator << (ostream& os,
		       const Trilinos_Util_MatrixGallery & G )
{

  bool verbose = (G.comm_->MyPID() == 0);

  if( verbose ) {
    
    os << " * Solving problem " << G.name_ << endl;
    os << " * Number of global elements : " << G.NumGlobalElements_ << endl;
    os << " * Type of Map : " << G.MapType_ << endl;
    os << " * Number of PDEs : " << G.NumPDEEqns_ << endl;

    // CRS stuff
    if( G.matrix_ != NULL ) {
      os << " * the matrix has been created " << endl;
      os << " * Matrix->OperatorDomainMap().NumGlobalElements() = "
	 << G.matrix_->OperatorDomainMap().NumGlobalElements() << endl;
    }
    if( G.ExactSolution_ != NULL ) 
      os << " * an exact solution (" << G.ExactSolutionType_
	 << ") has been created " << endl;
    if( G.rhs_ != NULL ) 
      os << " * the RHS has been created " << endl;
    // VBR stuff
    if( G.VbrMatrix_ != NULL ) {
      os << " * VBR matrix has been created " << endl;
      os << " * VBR Expanded as : " << G.MapType_ << endl;
      os << " * VbrMatrix->OperatorDomainMap().NumGlobalElements() = "
	 << G.VbrMatrix_->OperatorDomainMap().NumGlobalElements() << endl;
    }

    if( G.VbrExactSolution_ != NULL ) 
      os << " * a VBR exact solution has been created " << endl;
    if( G.VbrRhs_ != NULL ) 
      os << " * the VBR RHS has been created " << endl;
    
  }

  //  os << " * On proc " << G.comm_->MyPID() << " there are "
  //     << G.NumMyElements_ << " elements" << endl;

  return os;
  
}


