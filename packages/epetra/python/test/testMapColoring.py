#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#              PyTrilinos.Epetra: Python Interface to Epetra
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
    import Epetra
except ImportError:
    from PyTrilinos import Epetra
    print >>sys.stderr, "Using system-installed Epetra"

from   numpy  import *
import unittest

##########################################################################

class EpetraMapColoringTestCase(unittest.TestCase):
    "TestCase for Epetra_MapColorings"

    def setUp(self):
        self.comm        = Epetra.PyComm()
        self.myPID       = self.comm.MyPID()
        self.numProc     = self.comm.NumProc()
        self.mySize      = 4
        self.globalSize  = self.numProc * self.mySize
        self.map         = Epetra.Map(self.globalSize,0,self.comm)
        self.comm.Barrier()

    def tearDown(self):
        self.comm.Barrier()

    def testConstructor0(self):
        "Test Epetra.MapColoring BlockMap constructor"
        mc = Epetra.MapColoring(self.map)
        self.assertEqual(isinstance(mc,Epetra.MapColoring), True)

    def testConstructor1(self):
        "Test Epetra.MapColoring (BlockMap,int) constructor"
        mc = Epetra.MapColoring(self.map,8)
        self.assertEqual(isinstance(mc,Epetra.MapColoring), True)

    def testConstructor2(self):
        "Test Epetra.MapColoring (BlockMap,array) constructor"
        colors = range(self.map.NumMyElements())
        mc = Epetra.MapColoring(self.map,colors)
        self.assertEqual(isinstance(mc,Epetra.MapColoring), True)

    def testConstructor3(self):
        "Test Epetra.MapColoring (BlockMap,array,int) constructor"
        colors = range(self.map.NumMyElements())
        mc = Epetra.MapColoring(self.map,colors,31)
        self.assertEqual(isinstance(mc,Epetra.MapColoring), True)

    def testConstructor4(self):
        "Test Epetra.MapColoring copy constructor"
        mc1 = Epetra.MapColoring(self.map)
        mc2 = Epetra.MapColoring(mc1)
        self.assertEqual(isinstance(mc2,Epetra.MapColoring), True)

    def testSetGet(self):
        "Test Epetra.MapColoring set/get methods"
        mc = Epetra.MapColoring(self.map,0)
        for i in range(self.map.NumMyElements()):
            self.assertEqual(mc[i], 0)
        for i in range(self.map.NumMyElements()):
            mc[i] = i
            self.assertEqual(mc[i], i)

    def testNumColors(self):
        "Test Epetra.MapColoring NumColors method"
        mc = Epetra.MapColoring(self.map)
        self.assertEqual(mc.NumColors(), 1)
        n = self.map.NumMyElements()
        for i in range(n):
            mc[i] = i
        self.assertEqual(mc.NumColors(), n)

    def testMaxNumColors(self):
        "Test Epetra.MapColoring MaxNumColors method"
        mc = Epetra.MapColoring(self.map)
        self.assertEqual(mc.MaxNumColors(), 1)
        n = self.map.NumMyElements()
        for i in range(n):
            mc[i] = i
        self.assertEqual(mc.MaxNumColors(), n)

    def testListOfColors(self):
        "Test Epetra.MapColoring ListOfColors method"
        colors = range(self.mySize)
        mc     = Epetra.MapColoring(self.map,colors)
        result = mc.ListOfColors()
        self.assertEqual(len(result), self.mySize);
        for i in range(self.mySize):
            self.assertEqual(result[i], colors[i])

    def testDefaultColor(self):
        "Test Epetra.MapColoring DefaultColor method"
        mc = Epetra.MapColoring(self.map)
        self.assertEqual(mc.DefaultColor(), 0)
        dc = 1
        mc = Epetra.MapColoring(self.map,dc)
        self.assertEqual(mc.DefaultColor(), dc)

    def testNumElementsWithColor(self):
        "Test Epetra.MapColoring NumElementsWithColor method"
        colors = range(self.mySize)
        mc     = Epetra.MapColoring(self.map,colors)
        for i in range(self.mySize):
            self.assertEqual(mc.NumElementsWithColor(i), 1)

    def testColorLIDList(self):
        "Test Epetra.MapColoring ColorLIDList method"
        colors = range(self.mySize)
        mc     = Epetra.MapColoring(self.map,colors)
        for i in range(self.mySize):
            result = mc.ColorLIDList(i)
            self.assertEqual(len(result), 1);
            self.assertEqual(result[0], i)

    def testElementColors(self):
        "Test Epetra.MapColoring ElementColors method"
        colors = range(self.mySize)
        mc     = Epetra.MapColoring(self.map,colors)
        result = mc.ElementColors()
        self.assertEqual(len(result), self.mySize);
        for i in range(self.mySize):
            self.assertEqual(result[i], i)

    def testGenerateBlockMap(self):
        "Test Epetra.MapColoring GenerateBlockMap method"
        colors = range(self.mySize)
        mc     = Epetra.MapColoring(self.map,colors)
        for color in colors:
            blockMap = mc.GenerateBlockMap(color)
            result   = blockMap.MyGlobalElements()
            self.assertEqual(len(result), 1)
            self.assertEqual(result[0],self.map.GID(color))

    def testGenerateMap(self):
        "Test Epetra.MapColoring GenerateMap method"
        colors = range(self.mySize)
        mc     = Epetra.MapColoring(self.map,colors)
        for color in colors:
            map    = mc.GenerateMap(color)
            result = map.MyGlobalElements()
            self.assertEqual(len(result), 1)
            self.assertEqual(result[0],self.map.GID(color))

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraMapColoringTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
          "\n**************************\nTesting Epetra.MapColoring\n**************************\n"
    verbosity = 2 * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
