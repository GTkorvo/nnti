#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#         PyTrilinos.Amesos: Python Interface to Amesos
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

# Imports.  Users importing an installed version of PyTrilinos should use the
# "from PyTrilinos import ..." syntax.  Here, the setpath module adds the build
# directory, including "PyTrilinos", to the front of the search path.  We thus
# use "import ..." for Trilinos modules.  This prevents us from accidentally
# picking up a system-installed version and ensures that we are testing the
# build module.
import sys

try:
    import setpath
    import Epetra, Amesos
except ImportError:
    from PyTrilinos import Epetra, Amesos
    print >>sys.stderr, "Using system-installed Epetra, Amesos"

import unittest

####################################################################

# Swahili support
haveSwahili = "Newp_Jambo" in dir(Amesos)

####################################################################

class AmesosTestCase(unittest.TestCase):
    "TestCase class for Amesos module"

    def setUp(self):
        self.comm = Epetra.PyComm()
        self.comm.Barrier()

    def tearDown(self):
        # This will help tame the printing
        self.comm.Barrier()

    def testVersion(self):
        "Test Amesos Amesos_Version function"
        front   = "Amesos Version "
        version = Amesos.Amesos_Version()
        self.assertEquals(version[:len(front)], front)

    def testHelloConstructor(self):
        "Test Amesos Hello constructor"
        hello = Amesos.Newp_Hello(self.comm)
        # No assert here; just executing without an exception is enough

    def testHelloPrint(self):
        "Test Amesos Hello print operator"
        hello = Amesos.Newp_Hello(self.comm)
        if self.comm.MyPID() == 0:
            result = "This will print out one line for each of the " + \
                     str(self.comm.NumProc()) + " processes \n\n"
        else:
            result = ""
        result += "Hello.  I am process %d" % self.comm.MyPID()
        string = hello.__str__()
        self.assertEquals(string, result)

    if haveSwahili:
        def testJamboConstructor(self):
            "Test Amesos Jambo constructor"
            jambo = Amesos.Newp_Jambo(self.comm)
            # No assert here; just executing without an exception is enough

        def testJamboPrint(self):
            "Test Amesos Jambo print operator"
            jambo = Amesos.Newp_Jambo(self.comm)
            if self.comm.MyPID() == 0:
                result = "This will print out one line for each of the " + \
                         str(self.comm.NumProc()) + " processes \n\n"
            else:
                result = ""
            result += "Jambo.  I am process %d" % self.comm.MyPID()
            string = jambo.__str__()
            self.assertEquals(string, result)

####################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(AmesosTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
       "\n*******************\nTesting Amesos\n*******************\n"
    verbosity = 2 * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))[0]
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
