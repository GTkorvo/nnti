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

import unittest
from   Numeric    import *

##########################################################################

class EpetraCrsMatrixTestCase(unittest.TestCase):
    "TestCase class for Epetra CrsMatrix objects"

    def setUp(self):
        self.comm      = Epetra.PyComm()
        self.myPID     = self.comm.MyPID()
        self.numProc   = self.comm.NumProc()
        self.mySize    = 11
        self.size      = self.mySize * self.numProc
        self.indexBase = 0
        self.rowMap    = Epetra.Map(self.size, self.indexBase, self.comm)
        mge            = list(self.rowMap.MyGlobalElements())
        if self.indexBase not in mge: mge.append(mge[ 0]-1)
        if self.size-1    not in mge: mge.append(mge[-1]+1)
        mge.sort()
        self.colMap    = Epetra.Map(-1, mge, self.indexBase, self.comm)
        self.nipr      = ones(self.mySize)
        self.nipr[1:-1] += 1
        self.nipr[2:-2] += 1

    def tearDown(self):
        self.comm.Barrier()

    def fillMatrixGlobal(self,matrix,complete=True):
        n = self.size
        for lrid in range(matrix.NumMyRows()):
            grid = matrix.GRID(lrid)
            if grid == 0:
                values  = [2,-1]
                indices = [0, 1]
            elif grid == n-1:
                values  = [ -1,  2]
                indices = [n-2,n-1]
            else:
                values  = [    -1,   2,    -1]
                indices = [grid-1,grid,grid+1]
            result = matrix.InsertGlobalValues(grid,values,indices)
        if complete: matrix.FillComplete()

    def fillMatrixLocal(self,matrix,complete=True):
        n = self.size
        o = 0
        if self.myPID > 0: o = 1
        for lrid in range(matrix.NumMyRows()):
            grid = matrix.GRID(lrid)
            if grid == 0:
                values  = [2,-1]
                indices = [0, 1]
            elif grid == n-1:
                values  = [      -1,     2]
                indices = [lrid+o-1,lrid+o]
            else:
                values  = [      -1,     2,      -1]
                indices = [lrid+o-1,lrid+o,lrid+o+1]
            matrix.InsertMyValues(lrid,values,indices)
        if complete: matrix.FillComplete()

    def fillMatrixSetitem(self,matrix,complete=True):
        n = self.size
        for lrid in range(matrix.NumMyRows()):
            grid = matrix.GRID(lrid)
            if grid > 0:
                matrix[grid,grid-1] = -1
            matrix[grid,grid] = 2
            if grid < n-1:
                matrix[grid,grid+1] = -1
        if complete: matrix.FillComplete()

    def createGraph(self):
        graph = Epetra.CrsGraph(Epetra.Copy, self.rowMap, 3)
        n = self.size
        for lrid in range(graph.NumMyRows()):
            grid = graph.GRID(lrid)
            if   grid == 0  : indices = [0,1]
            elif grid == n-1: indices = [n-2,n-1]
            else            : indices = [grid-1,grid,grid+1]
            result = graph.InsertGlobalIndices(grid,indices)
        graph.FillComplete()
        return graph

    def createIdentity(self):
        ident = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 1, True)
        for lrid in range(ident.NumMyRows()):
            grid = ident.GRID(lrid)
            ident.InsertGlobalValues(grid,1.0,grid)
        ident.FillComplete()
        return ident

    def testConstructor01(self):
        "Test Epetra.CrsMatrix constructor, no colMap w/fixed # of indices/row"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.NumMyRows(), self.mySize   )
        self.assertEqual(crsm.IndexBase(), self.indexBase)

    def testConstructor02(self):
        "Test Epetra.CrsMatrix constructor, no colMap w/fixed # of indices/row, static profile"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3, True)
        self.assertEqual(crsm.NumMyRows(), self.mySize   )
        self.assertEqual(crsm.IndexBase(), self.indexBase)

    def testConstructor03(self):
        "Test Epetra.CrsMatrix constructor, no colMap w/specified # of indices/row"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.nipr)
        self.assertEqual(crsm.NumMyRows(), self.mySize   )
        self.assertEqual(crsm.IndexBase(), self.indexBase)

    def testConstructor04(self):
        "Test Epetra.CrsMatrix constructor, no colMap w/specified # of indices/row, static profile"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.nipr, True)
        self.assertEqual(crsm.NumMyRows(), self.mySize   )
        self.assertEqual(crsm.IndexBase(), self.indexBase)

    def testConstructor05(self):
        "Test Epetra.CrsMatrix constructor, w/colMap w/fixed # of indices/row"
        try:
            crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3)
            self.assertEqual(crsm.NumMyRows(), self.mySize   )
            self.assertEqual(crsm.IndexBase(), self.indexBase)
        except TypeError:
            # A TypeError is raised for the wrapper generated by older versions
            # of SWIG (1.3.28 and earlier).  This is expected, so just ignore it
            pass

    def testConstructor06(self):
        "Test Epetra.CrsMatrix constructor, w/colMap w/fixed # of indices/row, static profile"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, True)
        self.assertEqual(crsm.NumMyRows(), self.mySize   )
        self.assertEqual(crsm.IndexBase(), self.indexBase)

    def testConstructor07(self):
        "Test Epetra.CrsMatrix constructor, w/colMap w/specified # of indices/row"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, self.nipr)
        self.assertEqual(crsm.NumMyRows(), self.mySize   )
        self.assertEqual(crsm.IndexBase(), self.indexBase)

    def testConstructor08(self):
        "Test Epetra.CrsMatrix constructor, w/colMap w/specified # of indices/row, static profile"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, self.nipr, True)
        self.assertEqual(crsm.NumMyRows(), self.mySize   )
        self.assertEqual(crsm.IndexBase(), self.indexBase)

    def testConstructor09(self):
        "Test Epetra.CrsMatrix CrsGraph constructor"
        crsg = self.createGraph()
        crsm  = Epetra.CrsMatrix(Epetra.Copy,crsg)
        self.assertEqual(crsm.NumMyRows(), crsg.NumMyRows())
        self.assertEqual(crsm.IndexBase(), crsg.IndexBase())

    def testConstructor10(self):
        "Test Epetra.CrsMatrix copy constructor"
        crsm1 = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.nipr)
        self.fillMatrixGlobal(crsm1)
        crsm2 = Epetra.CrsMatrix(crsm1)
        self.assertEqual(crsm2.NumMyRows(), self.mySize   )
        self.assertEqual(crsm2.IndexBase(), self.indexBase)

    def testConstructor11(self):
        "Test Epetra.CrsMatrix constructor with too-short array"
        self.assertRaises(ValueError, Epetra.CrsMatrix, Epetra.Copy, self.rowMap,
                          self.nipr[1:-1])

    def testConstructor12(self):
        "Test Epetra.CrsMatrix constructor with too-long array"
        nipr = list(self.nipr)
        nipr.append(1)
        self.assertRaises(ValueError, Epetra.CrsMatrix, Epetra.Copy, self.rowMap,
                          nipr)

    def testPutScalar(self):
        "Test Epetra.CrsMatrix PutScalar method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        crsm.PutScalar(1.0)
        for lrid in range(crsm.NumMyRows()):
            grid = crsm.GRID(lrid)
            refValues = [1.0,1.0,1.0]
            if grid in (0,self.size-1): refValues.pop()
            (values,indices) = crsm.ExtractMyRowCopy(lrid)
            for i in range(len(values)):
                self.assertEqual(values[i],refValues[i])

    def testScale(self):
        "Test Epetra.CrsMatrix Scale method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixGlobal(crsm)
        crsm.Scale(10)
        for lrid in range(crsm.NumMyRows()):
            grid = crsm.GRID(lrid)
            refValues = [-10,20,-10]
            if grid == 0          : refValues.pop( 0)
            if grid == self.size-1: refValues.pop(-1)
            (values,indices) = crsm.ExtractMyRowCopy(lrid)
            for i in range(len(values)):
                self.assertEqual(values[i],refValues[i])

    def testInsertGlobalValues(self):
        "Test Epetra.CrsMatrix InsertGlobalValues method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)  # This calls crsm.InsertGlobalValues()
        n = self.size
        for lrid in range(crsm.NumMyRows()):
            grid = crsm.GRID(lrid)
            if grid in (0,n-1): numEntries = 2
            else:               numEntries = 3
            self.assertEqual(crsm.NumGlobalEntries(grid), numEntries)

    def testInsertGlobalValuesBad1(self):
        "Test Epetra.CrsMatrix InsertGlobalValues method for bad values"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        grid = crsm.GRID(0)
        self.assertRaises(TypeError, crsm.InsertGlobalValues, grid,
                          [0,"e","pi"], [grid,grid+1,grid+2])

    def testInsertGlobalValuesBad2(self):
        "Test Epetra.CrsMatrix InsertGlobalValues method for bad indices"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        grid = crsm.GRID(0)
        self.assertRaises(TypeError, crsm.InsertGlobalValues, grid,
                          [-1, 2, -1], [0,"e","pi"])

    def testReplaceGlobalValues(self):
        "Test Epetra.CrsMatrix ReplaceGlobalValues method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        grid      = crsm.GRID(1)
        indices   = [grid-1,grid,grid+1]
        newValues = [-10,20,-10]
        crsm.ReplaceGlobalValues(grid,newValues,indices)
        (values,indices) = crsm.ExtractGlobalRowCopy(grid)
        for i in range(len(values)):
            self.assertEqual(values[i], newValues[i])

    def testSumIntoGlobalValues(self):
        "Test Epetra.CrsMatrix SumIntoGlobalValues method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        grid = crsm.GRID(0)
        (values,indices) = crsm.ExtractGlobalRowCopy(grid)
        crsm.SumIntoGlobalValues(grid,values,indices)
        (newValues,indices) = crsm.ExtractGlobalRowCopy(grid)
        for i in range(len(values)):
            self.assertEqual(newValues[i], 2*values[i])

    def testInsertMyValues(self):
        "Test Epetra.CrsMatrix InsertMyValues method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixLocal(crsm)  # This calls crsm.InsertMyValues()
        n = self.size
        for lrid in range(crsm.NumMyRows()):
            grid = crsm.GRID(lrid)
            if grid in (0,n-1): numEntries = 2
            else              : numEntries = 3
            self.assertEqual(crsm.NumGlobalEntries(grid), numEntries)

    def testReplaceMyValues(self):
        "Test Epetra.CrsMatrix ReplaceMyValues method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixLocal(crsm,False)
        lrid      = 1
        if self.myPID > 0: o = 1
        else             : o = 0
        indices   = [lrid+o-1,lrid+o,lrid+o+1]
        newValues = [-10,20,-10]
        crsm.ReplaceMyValues(lrid,newValues,indices)
        (values,indices) = crsm.ExtractMyRowCopy(lrid)
        for i in range(len(values)):
            self.assertEqual(values[i], newValues[i])

    def testSumIntoMyValues(self):
        "Test Epetra.CrsMatrix SumIntoMyValues method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixLocal(crsm,False)
        lrid = 0
        (values,indices) = crsm.ExtractMyRowCopy(lrid)
        crsm.SumIntoMyValues(lrid,values,indices)
        (newValues,indices) = crsm.ExtractMyRowCopy(lrid)
        for i in range(len(values)):
            self.assertEqual(newValues[i], 2*values[i])

    def testReplaceDiagonalValues(self):
        "Test Epetra.CrsMatrix ReplaceDiagonalValues method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        ev = Epetra.Vector(self.rowMap)
        ev.PutScalar(10.0)
        crsm.ReplaceDiagonalValues(ev)
        for lrid in range(crsm.NumMyRows()):
            grid = crsm.GRID(lrid)
            (values,indices) = crsm.ExtractGlobalRowCopy(grid)
            refValues = [-1,10,-1]
            if grid == 0          : refValues.pop( 0)
            if grid == self.size-1: refValues.pop(-1)

    def testFillComplete1(self):
        "Test Epetra.CrsMatrix FillComplete method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.Filled(), False)
        result = crsm.FillComplete()
        self.assertEqual(result,        0   )
        self.assertEqual(crsm.Filled(), True)

    def testFillComplete2(self):
        "Test Epetra.CrsMatrix FillComplete method w/specified maps"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.Filled(), False)
        result = crsm.FillComplete(self.rowMap, self.colMap)
        self.assertEqual(result,        0   )
        self.assertEqual(crsm.Filled(), True)

    def testOptimizeStorage(self):
        "Test Epetra.CrsMatrix OptimizeStorage method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm,False)
        grid = crsm.GRID(0)
        numEntries = crsm.NumGlobalEntries(grid)
        crsm.InsertGlobalValues(grid,[2],[grid])  # Duplicate, thus non-optimal
        self.assertEqual(crsm.NumGlobalEntries(grid), numEntries+1)
        crsm.FillComplete()
        crsm.OptimizeStorage()
        self.assertEqual(crsm.NumGlobalEntries(grid), numEntries)

    def testMakeDataContiguous(self):
        "Test Epetra.CrsMatrix MakeDataContiguous method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm,False)
        grid = crsm.GRID(0)
        indices = range(grid+2,grid+5)
        values  = ones((3,),'d')
        crsm.InsertGlobalValues(grid,values,indices)
        self.assertEqual(crsm.IndicesAreContiguous(),False)
        crsm.FillComplete()
        crsm.MakeDataContiguous()
        self.assertEqual(crsm.IndicesAreContiguous(),True)

    def testExtractGlobalRowCopy(self):
        "Test Epetra.CrsMatrix ExtractGlobalRowCopy method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        n = crsm.GRID(self.mySize-1)
        (values,indices) = crsm.ExtractGlobalRowCopy(n)
        ncol = 3
        if n in (0,self.size-1): ncol -= 1
        self.assertEqual(len(values) , ncol)
        self.assertEqual(values[0]   , -1  )
        self.assertEqual(values[1]   , 2   )
        self.assertEqual(len(indices), ncol)
        self.assertEqual(indices[0]  , n-1 )
        self.assertEqual(indices[1]  , n   )

    def testExtractGlobalRowCopyBad(self):
        "Test Epetra.CrsMatrix ExtractGlobalRowCopy method, bad index"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        self.assertRaises(ValueError, crsm.ExtractGlobalRowCopy, self.size)

    def testExtractMyRowCopy(self):
        "Test Epetra.CrsMatrix ExtractMyRowCopy method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        n = self.mySize-1
        (values,indices) = crsm.ExtractMyRowCopy(n)
        ncol = 3
        if crsm.GRID(n) in (0,self.size-1): ncol -= 1
        self.assertEqual(len(values) , ncol)
        self.assertEqual(values[0]   , -1  )
        self.assertEqual(values[1]   , 2   )
        self.assertEqual(len(indices), ncol)
        self.assertEqual(indices[0]  , n-1 )
        self.assertEqual(indices[1]  , n   )

    def testExtractMyRowCopyBad(self):
        "Test Epetra.CrsMatrix ExtractMyRowCopy method, bad index"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        self.assertRaises(ValueError, crsm.ExtractMyRowCopy, self.mySize)

    def testExtractDiagonalCopy(self):
        "Test Epetra.CrsMatrix ExtractDiagonalCopy method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        ev = Epetra.Vector(self.rowMap)
        ev.PutScalar(-1.0)
        crsm.ExtractDiagonalCopy(ev)
        for i in range(len(ev)):
            self.assertEqual(ev[i],2.0)

    def testMultiplyVector(self):
        "Test Epetra CrsMatrix Multiply method for a Vector"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixGlobal(crsm)
        x = Epetra.Vector(self.rowMap)
        y = Epetra.Vector(self.rowMap)
        x.PutScalar(1.0)
        crsm.Multiply(False,x,y)   # y = A*x
        for lid in range(len(y)):
            gid = self.rowMap.GID(lid)
            result = 0.0
            if gid in (0, self.size-1): result = 1.0
            self.assertEqual(y[lid], result)

    def testMultiplyMultiVector(self):
        "Test Epetra CrsMatrix Multiply method for a MultiVector"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixGlobal(crsm)
        x = Epetra.MultiVector(self.rowMap,2)
        y = Epetra.MultiVector(self.rowMap,2)
        x.PutScalar(1.0)
        crsm.Multiply(False,x,y)   # y = A*x
        for lid in range(len(y)):
            gid = self.rowMap.GID(lid)
            result = 0.0
            if gid in (0, self.size-1): result = 1.0
            self.assertEqual(y[0,lid], result)
            self.assertEqual(y[1,lid], result)

    def testSolveVector(self):
        "Test Epetra CrsMatrix Solve method for a Vector"
        crsm = self.createIdentity()
        x = Epetra.Vector(self.rowMap)
        y = Epetra.Vector(self.rowMap)
        x.Random()
        result = crsm.Solve(True,False,False,x,y)
        self.assertEqual(result, 0)
        for i in range(len(x)):
            self.assertAlmostEqual(x[i], y[i])

    def testSolveMultiVector(self):
        "Test Epetra CrsMatrix Solve method for a MultiVector"
        crsm = self.createIdentity()
        x = Epetra.MultiVector(self.rowMap,2)
        y = Epetra.MultiVector(self.rowMap,2)
        x.Random()
        result = crsm.Solve(True,False,False,x,y)
        self.assertEqual(result, 0)
        for i in range(x.MyLength()):
            self.assertAlmostEqual(x[0,i], y[0,i])
            self.assertAlmostEqual(x[1,i], y[1,i])

    def testInvRowSums(self):
        "Test Epetra CrsMatrix InvRowSums method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        x = Epetra.Vector(self.rowMap)
        result = crsm.InvRowSums(x)
        self.assertEqual(result, 0)
        for lid in range(x.MyLength()):
            gid = self.rowMap.GID(lid)
            sum = 4.0
            if gid in (0,self.size-1): sum -= 1.0
            self.assertAlmostEqual(x[lid], 1.0/sum)

    def testInvRowMaxs(self):
        "Test Epetra CrsMatrix InvRowMaxs method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        x = Epetra.Vector(self.rowMap)
        result = crsm.InvRowMaxs(x)
        self.assertEqual(result, 0)
        for lid in range(x.MyLength()):
            self.assertAlmostEqual(x[lid], 1.0/2.0)

    def testLeftScale(self):
        "Test Epetra CrsMatrix LeftScale method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixGlobal(crsm)
        gid0 = self.rowMap.GID(0)
        x = Epetra.Vector(self.rowMap, arange(gid0+1, gid0+self.mySize+1))
        result = crsm.LeftScale(x)
        for lrid in range(crsm.NumMyRows()):
            grid = crsm.GRID(lrid)
            (values,indices) = crsm.ExtractGlobalRowCopy(grid)
            i = 0
            if grid > 0:
                self.assertEqual(values[i], -(grid+1))
                i += 1
            self.assertEqual(values[i], 2*(grid+1))
            i += 1
            if grid < self.size-1:
                self.assertEqual(values[i], -(grid+1))

    def testInvColSums(self):
        "Test Epetra CrsMatrix InvColSums method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        x = Epetra.Vector(self.rowMap)
        result = crsm.InvColSums(x)
        self.assertEqual(result, 0)
        for lid in range(x.MyLength()):
            gid = self.rowMap.GID(lid)
            sum = 4.0
            if gid in (0,self.size-1): sum -= 1.0
            self.assertAlmostEqual(x[lid], 1.0/sum)

    def testInvColMaxs(self):
        "Test Epetra CrsMatrix InvColMaxs method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        x = Epetra.Vector(self.rowMap)
        result = crsm.InvColMaxs(x)
        self.assertEqual(result, 0)
        for lid in range(x.MyLength()):
            gid = self.rowMap.GID(lid)
            self.assertAlmostEqual(x[lid], 1.0/2.0)

    def testRightScale(self):
        "Test Epetra CrsMatrix RightScale method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixGlobal(crsm)
        gid0 = self.rowMap.GID(0)
        x = Epetra.Vector(self.rowMap, arange(gid0+1, gid0+self.mySize+1))
        result = crsm.RightScale(x)
        for lrid in range(crsm.NumMyRows()):
            grid = crsm.GRID(lrid)
            (values,indices) = crsm.ExtractGlobalRowCopy(grid)
            i = 0
            if grid > 0:
                self.assertEqual(values[i], -grid)
                i += 1
            self.assertEqual(values[i], 2*(grid+1))
            i += 1
            if grid < self.size-1:
                self.assertEqual(values[i], -(grid+2))

    def testFilled(self):
        "Test Epetra.CrsMatrix Filled method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.Filled(), False)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.Filled(), True)

    def testStorageOptimized(self):
        "Test Epetra.CrsMatrix StorageOptimized method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.StorageOptimized(), False)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.StorageOptimized(), False)
        crsm.OptimizeStorage()
        self.assertEqual(crsm.StorageOptimized(), True)

    def testIndicesAreGlobal(self):
        "Test Epetra.CrsMatrix IndicesAreGlobal method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.IndicesAreGlobal(), False)
        self.fillMatrixGlobal(crsm,False)
        self.assertEqual(crsm.IndicesAreGlobal(), True)

    def testIndicesAreLocal(self):
        "Test Epetra.CrsMatrix IndicesAreLocal method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.assertEqual(crsm.IndicesAreLocal(), False)
        self.fillMatrixLocal(crsm)
        self.assertEqual(crsm.IndicesAreLocal(), True)

    def testIndicesAreContiguous(self):
        "Test Epetra.CrsMatrix IndicesAreContiguous method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        crsm.OptimizeStorage()
        self.assertEqual(crsm.IndicesAreContiguous(), True)

    def testLowerTriangular(self):
        "Test Epetra.CrsMatrix LowerTriangular method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.LowerTriangular(), True)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.LowerTriangular(), False)

    def testUpperTriangular(self):
        "Test Epetra.CrsMatrix UpperTriangular method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.UpperTriangular(), True)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.UpperTriangular(), False)

    def testNoDiagonal(self):
        "Test Epetra.CrsMatrix NoDiagonal method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.NoDiagonal(), True)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.NoDiagonal(), False)

    def testNormInf(self):
        "Test Epetra CrsMatrix NormInf method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.NormInf(), 4.0)

    def testNormOne(self):
        "Test Epetra CrsMatrix NormOne method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.NormOne(), 4.0)

    def testNumGlobalNonzeros(self):
        "Test Epetra.CrsMatrix NumGlobalNonzeros method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.NumGlobalNonzeros(), self.size*3-2)

    def testNumGlobalRows(self):
        "Test Epetra.CrsMatrix NumGlobalRows method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.NumGlobalRows(), self.size)

    def testNumGlobalCols(self):
        "Test Epetra.CrsMatrix NumGlobalCols method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.NumGlobalCols(), self.size)

    def testNumGlobalDiagonals(self):
        "Test Epetra.CrsMatrix NumGlobalDiagonals method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.NumGlobalDiagonals(), self.size)

    def testNumMyNonzeros(self):
        "Test Epetra.CrsMatrix NumMyNonzeros method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        numMyNonzeros = self.mySize * 3
        if crsm.MyGlobalRow(0          ): numMyNonzeros -= 1
        if crsm.MyGlobalRow(self.size-1): numMyNonzeros -= 1
        self.assertEqual(crsm.NumMyNonzeros(), numMyNonzeros)

    def testNumMyRows(self):
        "Test Epetra.CrsMatrix NumMyRows method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.NumMyRows(), self.mySize)

    def testNumMyCols(self):
        "Test Epetra.CrsMatrix NumMyCols method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        numMyCols = self.mySize
        if not crsm.MyGlobalRow(0          ): numMyCols += 1
        if not crsm.MyGlobalRow(self.size-1): numMyCols += 1
        self.assertEqual(crsm.NumMyCols(), numMyCols)

    def testNumMyDiagonals(self):
        "Test Epetra.CrsMatrix NumMyDiagonals method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.NumMyDiagonals(), self.mySize)

    def testNumGlobalEntries(self):
        "Test Epetra.CrsMatrix NumGlobalEntries method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        for lrid in range(crsm.NumMyRows()):
            grid = crsm.GRID(lrid)
            numEntries = 3
            if grid in (0, self.size-1): numEntries -= 1
            self.assertEqual(crsm.NumGlobalEntries(grid), numEntries)

    def testNumAllocatedGlobalEntries(self):
        "Test Epetra.CrsMatrix NumAllocatedGlobalEntries method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        for grid in self.rowMap.MyGlobalElements():
            self.assertEqual(crsm.NumAllocatedGlobalEntries(grid), 3)

    def testMaxNumEntries(self):
        "Test Epetra.CrsMatrix MaxNumEntries method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.MaxNumEntries(), 3)

    def testGlobalMaxNumEntries(self):
        "Test Epetra.CrsMatrix GlobalMaxNumEntries method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.GlobalMaxNumEntries(), 3)

    def testNumMyEntries(self):
        "Test Epetra.CrsMatrix NumMyEntries method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        for lrid in range(crsm.NumMyRows()):
            grid = crsm.GRID(lrid)
            numEntries = 3
            if grid in (0, self.size-1): numEntries -= 1
            self.assertEqual(crsm.NumMyEntries(lrid), numEntries)

    def testNumAllocatedMyEntries(self):
        "Test Epetra.CrsMatrix NumAllocatedMyEntries method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        for lrid in range(self.mySize):
            self.assertEqual(crsm.NumAllocatedMyEntries(lrid), 3)

    def testIndexBase(self):
        "Test Epetra.CrsMatrix IndexBase method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        self.assertEqual(crsm.IndexBase(), self.indexBase)

    def testStaticGraph(self):
        "Test Epetra.CrsMatrix StaticGraph method"
        crsm1 = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm1.StaticGraph(), False)
        crsg = self.createGraph()
        crsm2 = Epetra.CrsMatrix(Epetra.Copy, crsg)
        self.assertEqual(crsm2.StaticGraph(), True)

    def testGraph(self):
        "Test Epetra.CrsMatrix Graph method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        graph = crsm.Graph()
        self.assertEqual(crsm.Filled(),           graph.Filled()          )
        self.assertEqual(crsm.StorageOptimized(), graph.StorageOptimized())
        self.assertEqual(crsm.NumMyRows(),        graph.NumMyRows()       )
        self.assertEqual(crsm.NumMyCols(),        graph.NumMyCols()       )
        self.assertEqual(crsm.NumGlobalRows(),    graph.NumGlobalRows()   )
        self.assertEqual(crsm.NumGlobalCols(),    graph.NumGlobalCols()   )

    def testRowMap(self):
        "Test Epetra.CrsMatrix RowMap method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.RowMap().SameAs(self.rowMap), True)

    def testHaveColMap(self):
        "Test Epetra.CrsMatrix HaveColMap method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.HaveColMap(), False)
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.assertEqual(crsm.HaveColMap(), True)

    def testColMap(self):
        "Test Epetra.CrsMatrix ColMap method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.ColMap().SameAs(self.rowMap), True)
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.assertEqual(crsm.ColMap().SameAs(self.colMap), True)

    def testDomainMap(self):
        "Test Epetra.CrsMatrix DomainMap method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        domainMap = crsm.DomainMap()
        self.failUnless(self.rowMap.SameAs(domainMap))

    def testRangeMap(self):
        "Test Epetra.CrsMatrix RangeMap method"
        crsm1 = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        rangeMap1 = crsm1.RangeMap()
        self.failUnless(self.rowMap.SameAs(rangeMap1))
        crsm2 = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        rangeMap2 = crsm2.RangeMap()
        self.failUnless(self.rowMap.SameAs(rangeMap2))

    def testImporter(self):
        "Test Epetra.CrsMatrix Importer method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        importer = crsm.Importer()
        if self.numProc == 1: importerType = type(None)
        else                : importerType = Epetra.Import
        self.assertEqual(isinstance(importer, importerType), True)

    def testExporter(self):
        "Test Epetra.CrsMatrix Exporter method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        exporter = crsm.Exporter()
        self.assertEqual(isinstance(exporter, type(None)), True)

    def testComm(self):
        "Test Epetra.CrsMatrix Comm method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        comm = crsm.Comm()
        self.assertEqual(comm.NumProc(), self.numProc)
        self.assertEqual(comm.MyPID()  , self.myPID  )

    def testLRID(self):
        "Test Epetra.CrsMatrix LRID method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        for grid in range(self.size):
            if   grid <  self.mySize* self.myPID   : lrid = -1
            elif grid >= self.mySize*(self.myPID+1): lrid = -1
            else: lrid = grid - self.mySize*self.myPID
            self.assertEqual(crsm.LRID(grid), lrid)

    def testGRID(self):
        "Test Epetra.CrsMatrix GRID method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        for lrid in range(self.mySize):
            grid = lrid + self.mySize*self.myPID
            self.assertEqual(crsm.GRID(lrid), grid)
        self.assertEqual(crsm.GRID(self.mySize), -1)

    def testLCID(self):
        "Test Epetra.CrsMatrix LCID method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        start = self.mySize* self.myPID
        end   = self.mySize*(self.myPID+1)-1
        ghost1 = self.mySize
        if self.myPID == 0: ghost1 -= 1
        ghost2 = ghost1 + 1
        for gcid in range(self.size):
            if   gcid <  start-1: lcid = -1
            elif gcid == start-1: lcid = ghost1
            elif gcid == end+1  : lcid = ghost2
            elif gcid >  end+1  : lcid = -1
            else: lcid = gcid - start
            self.assertEqual(crsm.LCID(gcid), lcid)

    def testGCID(self):
        "Test Epetra.CrsMatrix GCID method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        size = self.mySize+2
        if self.myPID == 0             : size -= 1
        if self.myPID == self.numProc-1: size -= 1
        ghost = 0
        for lcid in range(size):
            if lcid == self.mySize:
                if self.myPID == 0: gcid = lcid
                else              : gcid = self.mySize*self.myPID-1
            elif lcid == self.mySize+1: gcid = lcid + self.mySize*self.myPID-1
            else                      : gcid = lcid + self.mySize*self.myPID
            self.assertEqual(crsm.GCID(lcid), gcid)
        self.assertEqual(crsm.GCID(self.mySize+2), -1)

    def testMyGRID(self):
        "Test Epetra.CrsMatrix MyGRID method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        for grid in range(self.size):
            if   grid <  self.mySize* self.myPID   : mine = False
            elif grid >= self.mySize*(self.myPID+1): mine = False
            else                                   : mine = True
            self.assertEqual(crsm.MyGRID(grid), mine)

    def testMyLRID(self):
        "Test Epetra.CrsMatrix MyLRID method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        for lrid in range(self.mySize):
            self.assertEqual(crsm.MyLRID(lrid), True)
        self.assertEqual(crsm.MyLRID(self.mySize), False)

    def testMyGCID(self):
        "Test Epetra.CrsMatrix MyGCID method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        start = self.mySize*self.myPID-1
        end = self.mySize*(self.myPID+1)
        for gcid in range(self.size):
            if   gcid < start: mine = False
            elif gcid > end  : mine = False
            else             : mine = True
            self.assertEqual(crsm.MyGCID(gcid), mine)

    def testMyLCID(self):
        "Test Epetra.CrsMatrix MyLCID method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixGlobal(crsm)
        size = self.mySize
        if self.myPID > 0             : size += 1
        if self.myPID < self.numProc-1: size += 1
        for lcid in range(size):
            self.assertEqual(crsm.MyLCID(lcid), True)
        self.assertEqual(crsm.MyLCID(size+1), False)

    def testMyGlobalRow(self):
        "Test Epetra.CrsMatrix MyGlobalRow method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        myRows = self.rowMap.MyGlobalElements()
        for grid in range(self.size):
            self.assertEqual(crsm.MyGlobalRow(grid), grid in myRows)

    def testLabel(self):
        "Test Epetra.CrsMatrix Label method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.Label(), "Epetra::CrsMatrix")

    def testSetUseTranspose(self):
        "Test Epetra.CrsMatrix SetUseTranspose and UseTranspose methods"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.UseTranspose(), False)
        crsm.SetUseTranspose(True)
        self.assertEqual(crsm.UseTranspose(), True)
        crsm.SetUseTranspose(False)
        self.assertEqual(crsm.UseTranspose(), False)

    def testApply(self):
        "Test Epetra CrsMatrix Apply method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixGlobal(crsm)
        x  = Epetra.Vector(self.rowMap)
        y1 = Epetra.Vector(self.rowMap)
        y2 = Epetra.Vector(self.rowMap)
        x.Random()
        crsm.Multiply(False,x,y1)   # y1 = A*x
        crsm.Apply(x,y2)            # y2 = A*x
        for i in range(x.MyLength()):
            self.assertEqual(y1[i], y2[i])

    def testApplyInverse(self):
        "Test Epetra CrsMatrix ApplyInverse method for diagonal matrix"
        crsm = self.createIdentity()
        x = Epetra.Vector(self.rowMap)
        y = Epetra.Vector(self.rowMap)
        x.Random()
        result = crsm.ApplyInverse(x,y)
        self.assertEqual(result, 0)
        for lid in range(x.MyLength()):
            self.assertEqual(x[lid], y[lid])

    def testApplyInverseBad(self):
        "Test Epetra CrsMatrix ApplyInverse method for non-triangular matrix"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        crsm.SetTracebackMode(0)
        self.fillMatrixGlobal(crsm)
        x = Epetra.Vector(self.rowMap)
        y = Epetra.Vector(self.rowMap)
        x.Random()
        result = crsm.ApplyInverse(x,y)
        self.assertEqual(result, -3)

    def testHasNormInf(self):
        "Test Epetra.CrsMatrix HasNormInf methods"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.assertEqual(crsm.HasNormInf(), True)

    def testOperatorDomainMap(self):
        "Test Epetra.CrsMatrix OperatorDomainMap method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        domainMap = crsm.OperatorDomainMap()
        self.failUnless(self.rowMap.SameAs(domainMap))

    def testOperatorRangeMap(self):
        "Test Epetra.CrsMatrix OperatorRangeMap method"
        crsm1 = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        rangeMap1 = crsm1.OperatorRangeMap()
        self.failUnless(self.rowMap.SameAs(rangeMap1))
        crsm2 = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        rangeMap2 = crsm2.OperatorRangeMap()
        self.failUnless(self.rowMap.SameAs(rangeMap2))

    def testSetitem(self):
        "Test Epetra CrsMatrix __setitem__ method"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        self.fillMatrixSetitem(crsm)
        for lrid in range(crsm.NumMyRows()):
            grid = crsm.GRID(lrid)
            if grid > 0: self.assertEqual(crsm[grid,grid-1], -1)
            self.assertEqual(crsm[grid,grid], 2)
            if grid < self.size-1: self.assertEqual(crsm[grid,grid+1], -1)

    def testSetitemBad(self):
        "Test Epetra CrsMatrix __setitem__ method with bad argument"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 3)
        grid = crsm.GRID(0)
        self.assertRaises(IndexError, crsm.__setitem__, grid, 1.0)

    def testGetitem1Global(self):
        "Test Epetra CrsMatrix __getitem__ method for 1 index, non-Filled"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixGlobal(crsm,False)
        grid0 = crsm.GRID(0)
        if self.numProc > 1:
            self.assertRaises(IndexError, crsm.__getitem__, grid0)
        else:
            for grid in range(grid0, grid0+crsm.NumMyRows()):
                rowData = crsm[grid]
                self.assertEqual(len(rowData), crsm.NumMyCols())
                for lcid in range(crsm.NumMyCols()):
                    gcid = crsm.GCID(lcid)
                    if gcid == grid-1:
                        self.assertEqual(rowData[lcid], -1.0)
                    elif gcid == grid:
                        self.assertEqual(rowData[lcid],  2.0)
                    elif gcid == grid+1:
                        self.assertEqual(rowData[lcid], -1.0)
                    else:
                        self.assertEqual(rowData[lcid],  0.0)

    def testGetitem1GlobalBad(self):
        "Test Epetra CrsMatrix __getitem__ method for 1 index, non-Filled, bad row ID"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixGlobal(crsm,False)
        grid0 = crsm.GRID(0)
        self.assertRaises(IndexError, crsm.__getitem__, grid0-1)

    def testGetitem1Local(self):
        "Test Epetra CrsMatrix __getitem__ method for 1 index, Filled"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixLocal(crsm)
        grid0 = crsm.GRID(0)
        for grid in range(grid0, grid0+crsm.NumMyRows()):
            rowData = crsm[grid]
            self.assertEqual(len(rowData), crsm.NumMyCols())
            for lcid in range(crsm.NumMyCols()):
                gcid = crsm.GCID(lcid)
                if gcid == grid-1:
                    self.assertEqual(rowData[lcid], -1.0)
                elif gcid == grid:
                    self.assertEqual(rowData[lcid],  2.0)
                elif gcid == grid+1:
                    self.assertEqual(rowData[lcid], -1.0)
                else:
                    self.assertEqual(rowData[lcid],  0.0)

    def testGetitem1LocalBad(self):
        "Test Epetra CrsMatrix __getitem__ method for 1 index, Filled, bad row ID"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixLocal(crsm)
        grid0 = crsm.GRID(0)
        self.assertRaises(IndexError, crsm.__getitem__, grid0-1)

    def testGetitem2Global(self):
        "Test Epetra CrsMatrix __getitem__ method for 2 indices, non-Filled"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixGlobal(crsm,False)
        grid0 = crsm.GRID(0)
        gcid0 = crsm.GCID(0)
        for grid in range(grid0, grid0+crsm.NumMyRows()):
            for gcid in range(gcid0, gcid0+crsm.NumMyCols()):
                if gcid == grid-1:
                    self.assertEqual(crsm[grid,gcid], -1.0)
                elif gcid == grid:
                    self.assertEqual(crsm[grid,gcid],  2.0)
                elif gcid == grid+1:
                    self.assertEqual(crsm[grid,gcid], -1.0)
                else:
                    self.assertEqual(crsm[grid,gcid],  0.0)

    def testGetitem2GlobalBad1(self):
        "Test Epetra CrsMatrix __getitem__ method for 2 indices, non-Filled, bad row ID"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixGlobal(crsm,False)
        grid0 = crsm.GRID(0)
        gcid0 = crsm.GCID(0)
        self.assertRaises(IndexError, crsm.__getitem__, (grid0-1, gcid0))

    def testGetitem2GlobalBad2(self):
        "Test Epetra CrsMatrix __getitem__ method for 2 indices, non-Filled, bad column ID"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixGlobal(crsm,False)
        grid0 = crsm.GRID(0)
        gcid0 = crsm.GCID(0)
        self.assertRaises(IndexError, crsm.__getitem__, (grid0, gcid0-1))

    def testGetitem2Local(self):
        "Test Epetra CrsMatrix __getitem__ method for 2 indices, Filled"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixLocal(crsm)
        grid0 = crsm.GRID(0)
        gcid0 = crsm.GCID(0)
        for grid in range(grid0, grid0+crsm.NumMyRows()):
            for gcid in range(gcid0, gcid0+crsm.NumMyCols()):
                if gcid == grid-1:
                    self.assertEqual(crsm[grid,gcid], -1.0)
                elif gcid == grid:
                    self.assertEqual(crsm[grid,gcid],  2.0)
                elif gcid == grid+1:
                    self.assertEqual(crsm[grid,gcid], -1.0)
                else:
                    self.assertEqual(crsm[grid,gcid],  0.0)

    def testGetitem2LocalBad1(self):
        "Test Epetra CrsMatrix __getitem__ method for 2 indices, Filled, bad row ID"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixLocal(crsm)
        grid0 = crsm.GRID(0)
        gcid0 = crsm.GCID(0)
        self.assertRaises(IndexError, crsm.__getitem__, (grid0-1, gcid0))

    def testGetitem2LocalBad2(self):
        "Test Epetra CrsMatrix __getitem__ method for 2 indices, Filled, bad column ID"
        crsm = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, self.colMap, 3, False)
        self.fillMatrixLocal(crsm)
        grid0 = crsm.GRID(0)
        gcid0 = crsm.GCID(0)
        self.assertRaises(IndexError, crsm.__getitem__, (grid0, gcid0-1))

##########################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(EpetraCrsMatrixTestCase))

    # Create a communicator
    comm    = Epetra.PyComm()
    iAmRoot = comm.MyPID() == 0

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
       "\n***********************\nTesting Epetra.CrsMatrix\n***********************\n"
    verbosity = 2 * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))[0]
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
