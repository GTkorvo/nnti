#! /usr/bin/env python

# @HEADER
# ************************************************************************
# 
#                PyTrilinos: Python Interface to Trilinos
#                   Copyright (2005) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Bill Spotz (wfspotz@sandia.gov) 
# 
# ************************************************************************
# @HEADER

from   optparse import *
import sys

parser = OptionParser()
parser.add_option("-t", "--testharness", action="store_true",
                  dest="testharness", default=False,
                  help="test local build modules; prevent loading system-installed modules")
parser.add_option("-v", "--verbosity", type="int", dest="verbosity", default=2,
                  help="set the verbosity level [default 2]")
options,args = parser.parse_args()
if options.testharness:
  import setpath
  import Epetra
  import Galeri
  import AztecOO
  import ML
else:
  try:
    import setpath
    import Epetra
    import Galeri
    import AztecOO
    import ML
  except ImportError:
    from PyTrilinos import Epetra
    from PyTrilinos import Galeri
    from PyTrilinos import AztecOO
    from PyTrilinos import ML
    print >>sys.stderr, "Using installed versions of Epetra, Galeri, AztecOO, ML"

def main(comm):

  # builds the linear system matrix and sets up starting solution and
  # right-hand side
  nx = 100
  ny = 100 * comm.NumProc()

  List = {
    "nx": nx,               # number of nodes in the X-direction
    "ny": ny,               # number of nodes in the Y-directioN
    "mx": 1,                # number of processors in the X-direction
    "my": comm.NumProc()    # number of processors in the Y-direction
  }

  Map = Galeri.CreateMap("Cartesian2D", comm, List)
  Matrix = Galeri.CreateCrsMatrix("Laplace2D", Map, List)

  LHS = Epetra.Vector(Map); LHS.Random()
  RHS = Epetra.Vector(Map); RHS.PutScalar(0.0)

  # sets up the parameters for ML using a python dictionary
  MLList = {
    "max levels"        : 3, 
    "output"            : 10,
    "smoother: type"    : "symmetric Gauss-Seidel",
    "aggregation: type" : "Uncoupled"
    }

  # creates the preconditioner and computes it
  Prec = ML.MultiLevelPreconditioner(Matrix, False)
  Prec.SetParameterList(MLList)
  Prec.ComputePreconditioner()

  # sets up the solver, specifies Prec as preconditioner, and
  # solves using CG.
  Solver = AztecOO.AztecOO(Matrix, LHS, RHS)
  Solver.SetPrecOperator(Prec)
  Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg);
  Solver.SetAztecOption(AztecOO.AZ_output, 16);
  err = Solver.Iterate(1550, 1e-5)

  return err

if __name__ == "__main__":
  comm = Epetra.PyComm()
  err = main(comm)
  errs = comm.SumAll(err)
  if errs == 0 and comm.MyPID() == 0: print "End Result: TEST PASSED"
  sys.exit(errs)
