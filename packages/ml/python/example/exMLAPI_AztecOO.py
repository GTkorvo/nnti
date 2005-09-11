#! /usr/bin/env python

# --------------------------------------------------------------------------- #
# Simple class that defines a multilevel preconditioner based on aggregation.
# The class requires in input the matrix (serial or distributed) defined as an
# ML.Operator, and the maximum number of levels. 
#
# NOTE: This code does not check for the actual number of levels; besides,
#       the smoother is always symmetric Gauss-Seidel. Simple changes can
#       be made to increase the flexibility of the code. If you want to 
#       define your own smoother, check example exMLAPI_Smoother.py
#
# \author Marzio Sala, SNL 9214
#
# \date Last updated on 03-Aug-05
# --------------------------------------------------------------------------- #

try:
  import setpath
except:
  from PyTrilinos import Epetra, ML, AztecOO
  print "Using installed versions of Epetra, ML, AztecOO"
else:
  import Epetra
  import ML
  import AztecOO

class MultiLevel(ML.BaseOperator):
  def Compute(self, A, MaxLevels):
    
    self.SetLabel("My MultiLevel Preconditioner")
    # Declares the lists that will contain the matrices (A), the prolongator
    # operators (P), the restriction operator (R), and the smoother or
    # coarse solver (S). The finest level matrix will be stored in A[0].
    A_array = [A]  
    P_array = []
    R_array = []
    S_array = []
    NullSpace = ML.MultiVector(A.GetRangeSpace())
    NullSpace.Update(1.0)

    self.MaxLevels_   = MaxLevels;
    self.A_           = A
    self.DomainSpace_ = A.GetDomainSpace()
    self.RangeSpace_  = A.GetRangeSpace()

    # Specify here the parameters for aggregation and smoothers
    List = {
      "aggregation: type": "Uncoupled",
      "aggregation: threshold": 0.05,
      "relaxation: sweeps": 2,
      "relaxation: damping factor": 0.67
    }

    # Builds the multilevel hierarchy
    for level in xrange(MaxLevels):
      A = A_array[level]
      NextNullSpace = ML.MultiVector()
      # Constructs the non-smoothed prolongator...
      Ptent = ML.GetPNonSmoothed(A, NullSpace, NextNullSpace, List)
      # ...and smooth it (without diagonal scaling, which is
      # easy to introduce if necessary)
      I = ML.GetIdentity(A.GetDomainSpace(), A.GetRangeSpace())
      lambda_max = 1.0 / ML.MaxEigAnorm(A)
      P = (I - A * lambda_max) * Ptent
      # Stores prolongator, restriction, and RAP product
      P = (I - A * lambda_max) * Ptent
      R = ML.GetTranspose(P)
      A_coarse = ML.GetRAP(R, A, P)
      # Defines the coarse solver or smoother (symmetric Gauss-Seidel)
      S = ML.InverseOperator()
      if level == MaxLevels - 1:
        Type = "Amesos"
      else:
        Type = "Jacobi"
      S.Reshape(A, Type, List)
    
      NullSpace = NextNullSpace
      A_array.append(A_coarse)
      P_array.append(P)
      R_array.append(R)
      S_array.append(S)

    self.A_array_ = A_array
    self.P_array_ = P_array
    self.R_array_ = R_array
    self.S_array_ = S_array

  def Apply(*args):
    self = args[0]
    RHS = args[1]
    LHS = args[2]
    LHS.Update(self.MultiLevelCycle(RHS, 0))
    return(0)
         
  def MultiLevelCycle(self, b_f, level):
    A = self.A_array_[level];
    P = self.P_array_[level];
    R = self.R_array_[level];
    S = self.S_array_[level];
    MaxLevels = self.MaxLevels_

    if level == MaxLevels - 1:
      return(S * b_f)

    # apply pre-smoother
    x_f = S * b_f
    # new residual
    r_f = b_f - A * x_f
    # restrict to coarse
    r_c = R * r_f
    # solve coarse problem
    z_c = self.MultiLevelCycle(r_c, level + 1)
    # prolongate back and add to solution
    x_f = x_f + P * z_c
    # apply post-smoother
    S.Apply(b_f, x_f)
  
    return(x_f)

  def GetOperatorDomainSpace(self):
    return(self.DomainSpace_)

  def GetOperatorRangeSpace(self):
    return(self.RangeSpace_)

  def __str__(self):
    return "MultiLevel"

  def __mul__(*args):
    self = args[0]
    rhs = args[1]
    res = ML.MultiVector(rhs.GetVectorSpace())
    self.Apply(rhs, res)
    return(res)

# -------------------------------------------------------------------------- #
# Main driver.
#
# The main difficulty of this example is that AztecOO requires Epetra
# objects (matrix, vectors, and preconditioner), while MLAPI is not
# an Epetra-derived class. However, it is easy (and computationally cheap)
# to wrap MLAPI object as Epetra object, and viceversa. In the following,
# you will have `_E' and `_M', meaning Epetra and MLAPI objects,
# respectively. Basically, you have to wrap the matrix and the
# preconditioner. The preconditioner class must be based on MLAPI only.
# -------------------------------------------------------------------------- #
#
def main():

  # Defines a communicator (serial or parallel, depending on how Trilinos
  # was configured), and creates a matrix corresponding to a 1D Laplacian.
  Comm = Epetra.PyComm()
  n = 1000
  Space = ML.Space(n)
  MyGlobalElements = Space.GetMyGlobalElements()
  
  Matrix = ML.PyMatrix(Space, Space)
  for i in MyGlobalElements:
    if i > 0:
      Matrix[i, i - 1] = -1.
    if i < n - 1:
      Matrix[i, i + 1] = -1.
    Matrix[i, i] = 2.0
  
  Matrix.FillComplete()
  
  MaxLevels = 3
  Prec = MultiLevel()
  Prec.Compute(Matrix, MaxLevels)
  
  # Both `Matrix' and `Prec' are MLAPI objects, that AztecOO cannot directly
  # understand. You must first extract a map for Matrix, then wrap the 
  # objects as Epetra objects
  Map = Matrix.GetMatrix().RowMatrixRowMap()
  MatrixWrap = ML.EpetraBaseOperator(Map, Matrix)
  # Defines a linear system, and solve it using AztecOO.
  LHS = Epetra.Vector(Map)
  RHS = Epetra.Vector(Map)
  RHS.Random()
  PrecWrap = ML.EpetraBaseOperator(Map, Prec)
 
  Solver = AztecOO.AztecOO(MatrixWrap, LHS, RHS)
  Solver.SetPrecOperator(PrecWrap)
  Solver.Iterate(10, 1e-5)

# This is a standard Python construct.  Put the code to be executed in a
# function [typically main()] and then use the following logic to call the
# function if the script has been called as an executable from the UNIX
# command
# line.  This also allows, for example, this file to be imported from a python
# debugger and main() called from there.
if __name__ == "__main__":
  main()
