
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */


#include "Epetra_SerialDenseVector.h"
//=============================================================================
Epetra_SerialDenseVector::Epetra_SerialDenseVector()
  : Epetra_SerialDenseMatrix()
{
	SetLabel("Epetra::SerialDenseVector");
}

//=============================================================================
Epetra_SerialDenseVector::Epetra_SerialDenseVector(int Length)
  : Epetra_SerialDenseMatrix(Length, 1)
{
	SetLabel("Epetra::SerialDenseVector");
}

//=============================================================================
Epetra_SerialDenseVector::Epetra_SerialDenseVector(Epetra_DataAccess CV, double *Values, int Length)
  : Epetra_SerialDenseMatrix(CV, Values, Length, Length, 1)
{
	SetLabel("Epetra::SerialDenseVector");	
}

//=============================================================================
Epetra_SerialDenseVector::Epetra_SerialDenseVector(const Epetra_SerialDenseVector& Source)
  : Epetra_SerialDenseMatrix(Source)
{}

//=============================================================================
Epetra_SerialDenseVector::~Epetra_SerialDenseVector()
{}

//=========================================================================
/*double& Epetra_SerialDenseVector::operator() (int Index)  {
  if (Index >= M_ || Index < 0)
		throw ReportError("Index = " +toString(Index) + " Out of Range 0 - " + toString(M_-1), -1); 
  return(A_[Index]); 
} 
//=========================================================================
const double& Epetra_SerialDenseVector::operator() (int Index) const  { 
  if (Index >= M_ || Index < 0)  
		throw ReportError("Index = " +toString(Index) + " Out of Range 0 - " + toString(M_-1), -1); 
   return(A_[Index]); 
}
//=========================================================================
const double& Epetra_SerialDenseVector::operator [] (int Index) const  { 
   return(A_[Index]); 
} 
//=========================================================================
double& Epetra_SerialDenseVector::operator [] (int Index)  { 
  if (Index >= M_ || Index < 0)  
		throw ReportError("Index = " +toString(Index) + " Out of Range 0 - " + toString(M_-1), -1); 
   return(A_[Index]); 
}*/

//=========================================================================
Epetra_SerialDenseVector& Epetra_SerialDenseVector::operator = (const Epetra_SerialDenseVector& Source) {
	Epetra_SerialDenseMatrix::operator=(Source); // call this->Epetra_SerialDenseMatrix::operator =
	return(*this);
}

//=========================================================================
int Epetra_SerialDenseVector::Random() {
	int errorcode = Epetra_SerialDenseMatrix::Random();
	return(errorcode);
}

//=========================================================================
void Epetra_SerialDenseVector::Print(ostream& os) const {
	if(CV_ == Copy)
		os << "Data access mode: Copy" << endl;
	else
		os << "Data access mode: View" << endl;
	if(A_Copied_)
		os << "A_Copied: yes" << endl;
	else
		os << "A_Copied: no" << endl;
	os << "Length(M): " << M_ << endl;
	if(M_ == 0)
		os << "(vector is empty, no values to display)";
	else
		for(int i = 0; i < M_; i++)
      os << (*this)(i) << " ";
	os << endl;
}
