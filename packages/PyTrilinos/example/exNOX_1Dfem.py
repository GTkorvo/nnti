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
# Questions? Contact Michael A. Heroux (maherou@sandia.gov)
#
# ************************************************************************
# @HEADER

#############################################################################
#
# This example script is intended to emulate the 1D FEM test problem from NOX
# found in nox/test/epetra/1Dfem/1Dfem.C.  It solves the 1D finite element
# approximation to
#
#     d2u
#     --- - k * u**2 = 0
#     dx2
#
# Subject to Dirichlet boundary conditions at x=0, x=1.
#
#############################################################################

# System imports
from   numpy    import *
from   optparse import *
import sys

parser = OptionParser()
parser.add_option("-t", "--testharness", action="store_true",
                  dest="testharness", default=False,
                  help="test local build modules; prevent loading system-installed modules")
parser.add_option("-v", "--verbosity", type="int", dest="verbosity", default=2,
                  help="set the verbosity level [default 2]")
parser.add_option("-n", "--numelem", type="int", dest="numelem", default=100,
                  help="set the number of elements [default 100]")
parser.add_option("-k", "--force_const", type="float", dest="k", default=1000.0,
                  help="set the nonlinear forcing constant [default 1000]")
parser.add_option("--bc0", type="float", dest="bc0", default=1.0,
                  help="set the Dirichlet boundary condition value at x=0 [default 1]")
parser.add_option("--bc1", type="float", dest="bc1", default=2.0,
                  help="set the Dirichlet boundary condition value at x=1 [default 2]")
parser.add_option("-p", "--plot", action="store_true", dest="plot",
                  default=False, help="Plot the resulting solution")
options,args = parser.parse_args()

# PyTrilinos imports.  If the -t or --testharness option is given, then only
# import the PyTrilinos modules from the current build directory structure.  If
# not, first try to load from the current build directories, but if that fails,
# try to import from a system-installed PyTrilinos.
if options.testharness:
    import setpath
    import Teuchos
    import Epetra
    import EpetraExt
    import NOX
else:
    try:
        import setpath
        import Teuchos
        import Epetra
        import EpetraExt
        import NOX
    except ImportError:
        from PyTrilinos import Teuchos
        from PyTrilinos import Epetra
        from PyTrilinos import EpetraExt
        from PyTrilinos import NOX
        print >>sys.stderr, "Using system-installed Teuchos, Epetra, EpetraExt, NOX"

# Plotting import
if options.plot:
    import pylab as mp

######################################################################

class Interface(NOX.Epetra.Interface.Required):

    def __init__(self, numGlobalElements, comm, bc0=1.0, bc1=2, xmin=0.0, xmax=1.0):
        # Initialize the base class
        NOX.Epetra.Interface.Required.__init__(self)
        # Interface initialization
        self.__numGlobalElements = numGlobalElements
        self.__comm              = comm
        self.__bc0               = bc0
        self.__bc1               = bc1
        self.__xmin              = xmin
        self.__xmax              = xmax
        self.__myPID             = comm.MyPID()
        self.__numProc           = comm.NumProc()
        self.__haveBC0           = (self.__myPID == 0)
        self.__haveBC1           = (self.__myPID == self.__numProc-1)
        self.__k                 = 1.0
        self.__stdMap            = Epetra.Map(numGlobalElements,0,comm)
        self.__numMyElements     = self.__stdMap.NumMyElements()
        newElements              = list(self.__stdMap.MyGlobalElements())
        if self.__myPID > 0:                newElements.insert(0,newElements[0]-1)
        if self.__myPID < self.__numProc-1: newElements.append(newElements[-1]+1)
        self.__overlapMap        = Epetra.Map(-1,newElements,0,comm)
        self.__importer          = Epetra.Import(self.__overlapMap, self.__stdMap)
        self.__initialSoln       = Epetra.Vector(self.__stdMap)
        self.__rhs               = Epetra.Vector(self.__stdMap)
        self.createGraph()
        x = array(self.__stdMap.MyGlobalElements()) / (numGlobalElements-1.0)
        self.__x                 = Epetra.Vector(self.__stdMap,x)
        self.__h                 = (self.__xmax - self.__xmin) / (len(self.__x)-1)
        self.initializeSoln()

    def computeF(self, u, F, flag):
        """
        This is the function that NOX calls back to in order to compute F(u).
        Arguments u and F are provided as Epetra.Vector objects, complete with
        numpy interface.
        """
        try:
            # Compute the residual on the interior
            F[1:-1] = (u[:-2] - 2*u[1:-1] + u[2:]) / (self.__h*self.__h) - \
                      self.__k * u[1:-1] * u[1:-1]
            # Compute the residual on the boundaries
            if self.__haveBC0: F[ 0] = u[ 0] - self.__bc0
            if self.__haveBC1: F[-1] = u[-1] - self.__bc1
            return True
        except Exception, e:
            print e
            return False

    def getSolution(self):
        return self.__initialSoln

    def getMesh(self):
        return self.__x

    def setPDEfactor(self, k):
        self.__k = k
        return True

    def createGraph(self):
        self.__graph = Epetra.CrsGraph(Epetra.Copy, self.__stdMap, 5)
        for lid in self.__stdMap.MyGlobalElements():
            gid = self.__stdMap.GID(lid)
            if gid in (0,self.__numGlobalElements-1):  # Boundaries
                self.__graph.InsertGlobalIndices(lid,[gid])
            else:                                      # Interior
                self.__graph.InsertGlobalIndices(lid,[gid-1,gid,gid+1])
        self.__graph.FillComplete()
        return True

    def initializeSoln(self):
        self.__initialSoln.PutScalar((self.__bc0 + self.__bc1)/2)
        if self.__haveBC0: self.__initialSoln[ 0] = self.__bc0
        if self.__haveBC1: self.__initialSoln[-1] = self.__bc1
        return True

    def getGraph(self):
        return self.__graph

######################################################################

# Main routine
def main():

    # Communicator
    comm    = Epetra.PyComm()
    myPID   = comm.MyPID()
    numProc = comm.NumProc()

    # Get the number of elements from the command line
    numGlobalElements = options.numelem + 1
    if numGlobalElements < numProc:
        msg = "numGlobalBlocks = %d cannot be < number of processors = %d" % \
              (numGlobalElements,numProc)
        msg += "\nTest failed!"
        raise RuntimeError, msg

    # Get the boundary conditions
    bc0 = options.bc0
    bc1 = options.bc1

    # Create the interface between NOX and the application.  This object is
    # derived from NOX.Epetra.Interface.
    interface = Interface(numGlobalElements, comm, bc0, bc1)

    # Get the solution vector from the problem
    soln    = interface.getSolution()
    noxSoln = NOX.Epetra.Vector(soln,NOX.Epetra.Vector.CreateView)

    # Set the nonlinear forcing constant
    interface.setPDEfactor(options.k)

    ### Begin Nonlinear Solve ###

    # Create the nonlinear solver parameters
    nlParams = {"Nonlinear Solver" : "Line Search Based",
                "Printing"         : {"MyPID"            : myPID,
                                      "Output Precision" : 3,
                                      "Output Processor" : 0
                                      },
                "Line Search"      : {"Method" : "Full Step"},
                "Direction"        : {"Method" : "Newton"},
                "Newton"           : {"Forcing Term Method" : "Constant"},
                "Linear Solver"    : {"Aztec Solver"    : "GMRES",
                                      "Max Iterations"  : 800,
                                      "Tolerance"       : 1e-4,
                                      "Preconditioner"  : "Ifpack",
                                      "Max Age Of Prec" : 5},
                #"Solver Options"   : {"Status Test Check Type" : NOX.StatusTest.Complete}
                }
    printParams = nlParams["Printing"]
    lsParams    = nlParams["Linear Solver"]
    outputInfo = NOX.Utils.Error + \
                 NOX.Utils.TestDetails
    if options.verbosity: outputInfo += NOX.Utils.Debug      + \
                                        NOX.Utils.Warning    + \
                                        NOX.Utils.Details    + \
                                        NOX.Utils.Parameters + \
                                        NOX.Utils.LinearSolverDetails
    if options.verbosity > 1: outputInfo += NOX.Utils.InnerIteration           + \
                                            NOX.Utils.OuterIterationStatusTest + \
                                            NOX.Utils.OuterIteration
    printParams["Output Information"] = outputInfo

    # Define the Jacobian linear system
    mf     = NOX.Epetra.MatrixFree(printParams,interface,noxSoln)
    fdc    = NOX.Epetra.FiniteDifferenceColoring(printParams, interface,
                                                 soln, interface.getGraph())
    linSys = NOX.Epetra.LinearSystemAztecOO(printParams, lsParams, mf, mf, fdc,
                                            fdc, soln)

    # Create the Group
    initialGuess = NOX.Epetra.Vector(soln, NOX.Epetra.Vector.CreateView)
    group        = NOX.Epetra.Group(printParams, interface, initialGuess, linSys)

    # Create the convergence tests
    absResid  = NOX.StatusTest.NormF(1.0e-8)
    relResid  = NOX.StatusTest.NormF(group,1.0e-2)
    update    = NOX.StatusTest.NormUpdate(1.0e-5)
    wrms      = NOX.StatusTest.NormWRMS(1.0e-2,1.0e-8)
    converged = NOX.StatusTest.Combo(NOX.StatusTest.Combo.AND)
    converged.addStatusTest(absResid)
    converged.addStatusTest(relResid)
    converged.addStatusTest(wrms)
    converged.addStatusTest(update)
    maxIters = NOX.StatusTest.MaxIters(20)
    fv       = NOX.StatusTest.FiniteValue()
    combo    = NOX.StatusTest.Combo(NOX.StatusTest.Combo.OR)
    combo.addStatusTest(fv       )
    combo.addStatusTest(converged)
    combo.addStatusTest(maxIters )

    # Create and execute the solver
    solver      = NOX.Solver.Manager(group, combo, nlParams)
    solveStatus = solver.solve()

    ### End Nonlinear Solver ###

    # Obtain the final solution
    finalGroup = solver.getSolutionGroup()
    finalSoln  = finalGroup.getX()

    # Output the parameter list
    if options.verbosity:
        if myPID == 0:
            print "Final Parameters\n****************"
            solver.getList()._print()
            print

    # Print solution
    #f = open("output." + str(myPID),"w")
    #finalSoln.Print(f)
    #close(f)

    # Tests
    status = 0
    if solveStatus != NOX.StatusTest.Converged:
        status += 1
        if myPID == 0: print "Nonlinear solver failed to converge"
    #if numProc == 0:
    #    if solver.getList().get("Direction").get("Newton").get("Linear Solver").get(
    #        "Output").get("Total Number of Linear Iterations") != 53:
    #        status += 1
    #if solver.getList().get("Output").get("Nonlinear Iterations") != 10:
    #    status += 1
    #if ppo.getNumRunPreIterate()  != 10: status += 1
    #if ppo.getNumRunPostIterate() != 10: status += 1
    #if ppo.getNumRunPreSolve()    !=  1: status += 1
    #if ppo.getNumRunPostSolve()   !=  1: status += 1

    # Summarize the test results
    if myPID == 0:
        if status == 0:
            print "End Result: TEST PASSED"
        else:
            print "End Result: TEST FAILED"

    # Plot the results
    if options.plot:
        if myPID == 0: print "Plotting results..."
        noxSoln.update(1.0,finalSoln,0.0)    # Copy the final solution into a NOX.Epetra.Vector
        f1 = noxSoln.getEpetraVector()       # Obtain the underlying Epetra.Epetra_Vector
        f2 = Epetra.Vector(Epetra.View,f1,0) # Convert to Epetra.Vector
        x  = interface.getMesh()
        mp.plot(x,f2)
        mp.title("exNOX_1Dfem Problem Solution, k = %d" % options.k)
        mp.xlabel("x")
        mp.ylabel("u(x)")
        mp.show()

    return status

######################################################################

if __name__ == "__main__":

    status = main()

    sys.exit(status)
