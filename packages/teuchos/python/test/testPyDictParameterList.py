#! /usr/bin/env python

# @HEADER
# ************************************************************************
#
#             PyTrilinos.Teuchos: Python Interface to Teuchos
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
from   optparse import *
import sys
import unittest

parser = OptionParser()
parser.add_option("-t", "--testharness", action="store_true",
                  dest="testharness", default=False,
                  help="test local build modules; prevent loading system-installed modules")
parser.add_option("-v", "--verbosity", type="int", dest="verbosity", default=2,
                  help="set the verbosity level [default 2]")
options,args = parser.parse_args()
if options.testharness:
    import setpath
    import Teuchos
else:
    try:
        import setpath
        import Teuchos
    except ImportError:
        from PyTrilinos import Teuchos
        print >>sys.stderr, "Using system-installed Teuchos"

####################################################################

class PyDictParameterListTestCase(unittest.TestCase):
    "TestCase class for Teuchos.PyDictParameterList"

    def setUp(self):
        self.plist = Teuchos.PyDictParameterList()
        self.name  = "Solver Params"

    def testConstructor0(self):
        "Test Teuchos.PyDictParameterList default constructor"
        self.assertEqual(isinstance(self.plist,Teuchos.PyDictParameterList), True)

    def testConstructor1(self):
        "Test Teuchos.PyDictParameterList string constructor"
        plist = Teuchos.PyDictParameterList(self.name)
        self.assertEqual(isinstance(plist,Teuchos.PyDictParameterList), True)
        self.assertEqual(plist.name(), self.name)

    def testConstructor2(self):
        "Test Teuchos.PyDictParameterList dictionary constructor"
        d = {"i" : 5,
             "f" : 3.14,
             "s" : "New Mexico"
             }
        plist = Teuchos.PyDictParameterList(d)
        for key in d.keys():
            self.assertEqual(plist.get(key), d[key])

    def testConstructor3(self):
        "Test Teuchos.PyDictParameterList nested dictionary constructor"
        d = {"i" : 5,
             "f" : 3.14,
             "s" : "New Mexico",
             "d" : {"a":1, "b":2}
             }
        plist = Teuchos.PyDictParameterList(d)
        for key in d.keys():
            if key == "d":
                sublist = plist.sublist(key)
                subdict = d["d"]
                for subKey in subdict.keys():
                    self.assertEqual(sublist.get(subKey), subdict[subKey])
            else:
                self.assertEqual(plist.get(key), d[key])

    def testConstructor4(self):
        "Test Teuchos.PyDictParameterList dictionary constructor, bad key"
        d = {"i" : 5,
             "f" : 3.14,
             "s" : "New Mexico",
             1   : "Should fail"
             }
        self.assertRaises(KeyError, Teuchos.PyDictParameterList, d)

    def testConstructor5(self):
        "Test Teuchos.PyDictParameterList dictionary constructor, bad value"
        d = {"i" : 5,
             "f" : 3.14,
             "s" : "New Mexico",
             "b" : (1,2,3)
             }
        self.assertRaises(TypeError, Teuchos.PyDictParameterList, d)

    def testConstructor6(self):
        "Test Teuchos.PyDictParameterList ParameterList constructor"
        plist = Teuchos.ParameterList()
        plist.set("unity",1)
        pyDictPlist = Teuchos.PyDictParameterList(plist)
        self.assertEqual(isinstance(pyDictPlist,Teuchos.PyDictParameterList), True)
        self.assertEqual(pyDictPlist.get("unity"),1)

    def testConstructor7(self):
        "Test Teuchos.PyDictParameterList copy constructor"
        plist_copy = Teuchos.PyDictParameterList(self.plist)
        self.assertEqual(isinstance(plist_copy,Teuchos.PyDictParameterList), True)

    def testSetName(self):
        "Test Teuchos.PyDictParameterList name and setName methods"
        self.assertEqual(self.plist.name(), "ANONYMOUS")
        self.plist.setName(self.name)
        self.assertEqual(self.plist.name(), self.name)

    def testSetParameters1(self):
        "Test Teuchos.PyDictParameterList setParameters method for a PyDictParameterList"
        intName    = "int parameter"
        intValue   = 8
        floatName  = "float parameter"
        floatValue = 3.14
        self.plist.set(intName,  intValue  )
        self.plist.set(floatName,floatValue)
        newList = Teuchos.PyDictParameterList()
        newList.setParameters(self.plist)
        self.assertEqual(newList.get(intName  ), intValue  )
        self.assertEqual(newList.get(floatName), floatValue)
        self.assertEqual(newList.isSynchronized(), True)

    def testSetParameters2(self):
        "Test Teuchos.PyDictParameterList setParameters method for a ParameterList"
        plist = Teuchos.ParameterList()
        plist.set("int parameter"  , 8)
        plist.set("float parameter", 3.14)
        self.plist.setParameters(plist)
        self.assertEqual(self.plist.isSynchronized(), True)

    def testSetParameters3(self):
        "Test Teuchos.PyDictParameterList setParameters method for a dictionary"
        d = {"int parameter"   : 8,
             "float parameter" : 3.14 }
        self.plist.setParameters(d)
        self.assertEqual(self.plist.isSynchronized(), True)

    def testSetInt(self):
        "Test Teuchos.PyDictParameterList set and get methods for an integer"
        name  = "int parameter"
        value = 12
        self.plist.set(name, value)
        self.assertEqual(self.plist.get(name), value)

    def testSetFloat(self):
        "Test Teuchos.PyDictParameterList set and get methods for a float"
        name  = "float parameter"
        value = 12.0
        self.plist.set(name, value)
        self.assertEqual(self.plist.get(name), value)
        self.assertEqual(self.plist.isSynchronized(), True)

    def testSetString(self):
        "Test Teuchos.PyDictParameterList set and get methods for a string"
        name  = "string parameter"
        value = "12"
        self.plist.set(name, value)
        self.assertEqual(self.plist.get(name), value)
        self.assertEqual(self.plist.isSynchronized(), True)

    def testSetDictionary(self):
        "Test Teuchos.PyDictParameterList set and get methods for a dictionary"
        name = "dict parameter"
        d    = {"a" : 1,
                "b" : 2,
                "c" : 3 }
        self.plist.set(name, d)
        sublist = self.plist.get(name)
        self.assertEqual(isinstance(sublist, Teuchos.ParameterList), True)
        for key in d.keys():
            self.assertEqual(sublist.get(key), d[key])
        self.assertEqual(self.plist.isSynchronized(), True)

    def testSetParameterList(self):
        "Test Teuchos.PyDictParameterList set and get methods for a ParameterList"
        name    = "sublist"
        sublist = Teuchos.ParameterList()
        self.plist.set(name, sublist)
        self.assertEqual(isinstance(self.plist.get(name),
                                    Teuchos.ParameterList), True)
        self.assertEqual(self.plist.isSublist(name), True)
        self.assertEqual(self.plist.isSynchronized(), True)

    def testSetPyDictParameterList(self):
        "Test Teuchos.PyDictParameterList set and get methods for a PyDictParameterList"
        name    = "sublist"
        sublist = Teuchos.PyDictParameterList()
        self.plist.set(name, sublist)
        self.assertEqual(isinstance(self.plist.get(name),
                                    Teuchos.ParameterList), True)
        self.assertEqual(self.plist.isSublist(name), True)
        self.assertEqual(self.plist.isSynchronized(), True)

    def testSetNone(self):
        "Test Teuchos.PyDictParameterList set method for None"
        self.assertRaises(ValueError, self.plist.set, "bad parameter", None)
        self.assertEqual(self.plist.isSynchronized(), True)

    def testSetBadType(self):
        "Test Teuchos.PyDictParameterList set method for an unsupported type"
        self.assertRaises(RuntimeError, self.plist.set, "bad parameter", (1,2,3))
        self.assertEqual(self.plist.isSynchronized(), True)

    def testGetDefault(self):
        "Test Teuchos.PyDictParameterList get method using default value"
        default = None
        self.assertEqual(self.plist.get("junk",None), None)

    def testSublistNew(self):
        "Test Teuchos.PyDictParameterList sublist method for new sublist"
        sublist = self.plist.sublist("new")
        self.assertEqual(isinstance(sublist, Teuchos.ParameterList), True)
        sublist = self.plist.get("new")
        self.assertEqual(isinstance(sublist, Teuchos.ParameterList), True)
        self.assertEqual(self.plist.isSynchronized(), True)

    def testSublistOld(self):
        "Test Teuchos.PyDictParameterList sublist method for existing sublist"
        sublist = self.plist.sublist("new")
        self.assertEqual(isinstance(sublist, Teuchos.ParameterList), True)
        sublist = self.plist.sublist("new")
        self.assertEqual(isinstance(sublist, Teuchos.ParameterList), True)
        self.assertEqual(self.plist.isSynchronized(), True)

    def testSublistBad(self):
        "Test Teuchos.PyDictParameterList sublist method for non-sublist"
        self.plist.set("new", 1)
        self.assertRaises(RuntimeError, self.plist.sublist, "new")

    def testIsParameterTrue(self):
        "Test Teuchos.PyDictParameterList isParameter method existing parameter"
        name = "string parameter"
        self.plist.set(name,"Hello")
        self.assertEqual(self.plist.isSynchronized(), True)
        self.assertEqual(self.plist.isParameter(name), True)

    def testIsParameterFalse(self):
        "Test Teuchos.PyDictParameterList isParameter method nonexisting parameter"
        name = "parameter"
        self.assertEqual(self.plist.isSynchronized(), True)
        self.assertEqual(self.plist.isParameter(name), False)

    def testIsSublistTrue(self):
        "Test Teuchos.PyDictParameterList isSublist method for existing sublist"
        name = "sublist"
        self.plist.sublist(name)
        self.assertEqual(self.plist.isSublist(name), True)

    def testIsSublistFalse1(self):
        "Test Teuchos.PyDictParameterList isSublist method for existing non-sublist parameter"
        name = "string parameter"
        self.plist.set(name,"Hello")
        self.assertEqual(self.plist.isSublist(name), False)

    def testIsSublistFalse2(self):
        "Test Teuchos.PyDictParameterList isSublist method for nonexisting parameter"
        name = "parameter"
        self.assertEqual(self.plist.isSublist(name), False)

    def testPrint0(self):
        "Test Teuchos.PyDictParameterList _print method for empty list"
        fName = "testPyDictParameterList.dat"
        f = open(fName, "w")
        self.plist._print(f)
        f.close()
        self.assertEqual(open(fName,"r").read(), "[empty list]\n")

    def testPrint1(self):
        "Test Teuchos.PyDictParameterList _print method for non-empty list"
        names  = ["max its","tolerance"]
        values = [100      , 1e-6      ]
        for i in range(len(names)):
            self.plist.set(names[i], values[i])
        fName = "testPyDictParameterList.dat"
        f = open(fName, "w")
        self.plist._print(f)
        f.close()
        lines = open(fName,"r").readlines()
        for i in range(len(lines)):
            self.assertEqual(lines[i], "%s = %g   [unused]\n" % (names[i], values[i]))

    def testPrint2(self):
        "Test Teuchos.PyDictParameterList _print method for non-empty list and indentation"
        names  = ["max its","tolerance"]
        values = [100      , 1e-6      ]
        for i in range(len(names)):
            self.plist.set(names[i], values[i])
        fName = "testPyDictParameterList.dat"
        f = open(fName, "w")
        self.plist._print(f,2)
        f.close()
        lines = open(fName,"r").readlines()
        for i in range(len(lines)):
            self.assertEqual(lines[i], "  %s = %g   [unused]\n" % (names[i], values[i]))

    def testPrint3(self):
        "Test Teuchos.PyDictParameterList _print method for non-empty list, indentation and types"
        names  = ["max its","tolerance"]
        values = [100      , 1e-6      ]
        types  = ["int"    ,"double"   ]
        for i in range(len(names)):
            self.plist.set(names[i], values[i])
        fName = "testPyDictParameterList.dat"
        f = open(fName, "w")
        self.plist._print(f,2,True)
        f.close()
        lines = open(fName,"r").readlines()
        for i in range(len(lines)):
            self.assertEqual(lines[i], "  %s : %s = %g   [unused]\n" %
                             (names[i], types[i], values[i]))

    def testUnused(self):
        "Test Teuchos.PyDictParameterList unused method"
        names  = ["s1"   , "s2"   , "s3"         ]
        values = ["Hello", "World", "Albuquerque"]
        for i in range(len(names)):
            self.plist.set(names[i], values[i])
        self.assertEqual(self.plist.isSynchronized(), True)
        fName = "testPyDictParameterList.dat"
        f = open(fName,"w")
        self.plist.unused(f)
        f.close()
        lines = open(fName,"r").readlines()
        for i in range(len(lines)):
            self.assertEqual(lines[i],
                             'WARNING: Parameter "%s" %s   [unused] is unused\n' %
                             (names[i], values[i]))

    def testCurrentParametersString(self):
        "Test Teuchos.PyDictParameterList currentParametersString method"
        names  = ["max its","tolerance"]
        values = [100      , 1e-6      ]
        types  = ["int"    ,"double"   ]
        result = "{"
        for i in range(len(names)):
            self.plist.set(names[i], values[i])
            result += '"%s":%s=%g,' % (names[i], types[i], values[i])
        result = result[:-1] + "}"
        self.assertEqual(self.plist.currentParametersString(), result)

    def testType(self):
        "Test Teuchos.PyDictParameterList type method"
        sublist = Teuchos.PyDictParameterList()
        names  = ["iParm", "fParm", "sParm", "lParm"              ]
        values = [2006   , 2.71828, "Hello", sublist              ]
        types  = [int    , float  , str    , Teuchos.ParameterList]
        for i in range(len(names)):
            self.plist.set(names[i],values[i])
        self.assertEqual(self.plist.isSynchronized(), True)
        for i in range(len(names)):
            self.assertEqual(self.plist.type(names[i]), types[i])

    def test__cmp__1(self):
        "Test Teuchos.PyDictParameterList __cmp__ method with dictionary"
        d1 = { "a" : 2 }
        d2 = { "a" : 2, "b" : 1 }
        plist = Teuchos.PyDictParameterList(d1)
        self.assertEqual(cmp(plist, d1), 0)
        self.assertEqual(cmp(plist, d2),-1)

    def test__cmp__2(self):
        "Test Teuchos.PyDictParameterList __cmp__ method with ParameterList"
        d = { "a" : 2, "b" : 1 }
        pdplist = Teuchos.PyDictParameterList(d)
        plist   = Teuchos.ParameterList(pdplist)
        self.assertEqual(cmp(pdplist, plist), 0)
        plist.set("c",3)
        self.assertEqual(cmp(pdplist, plist), 1)

    def test__contains__(self):
        "Test Teuchos.PyDictParameterList __contains__ method"
        d = { "a" : 2, "b" : 1 }
        plist = Teuchos.PyDictParameterList(d)
        self.assertEqual("b" in d, True )
        self.assertEqual("c" in d, False)

    def test__eq__1(self):
        "Test Teuchos.PyDictParameterList __eq__ method with dictionary"
        d1 = { "a" : 2 }
        d2 = { "a" : 2, "b" : 1 }
        plist = Teuchos.PyDictParameterList(d1)
        self.assertEqual(plist == d1, True )
        self.assertEqual(plist == d2, False)

    def test__eq__2(self):
        "Test Teuchos.PyDictParameterList __eq__ method with ParameterList"
        d = { "a" : 2, "b" : 1 }
        pdplist = Teuchos.PyDictParameterList(d)
        plist   = Teuchos.ParameterList(pdplist)
        self.assertEqual(pdplist == plist, True )
        plist.set("c",3)
        self.assertEqual(pdplist == plist, False)

    def test__getitem__(self):
        "Test Teuchos.PyDictParameterList __getitem__ method"
        d = {"i" : 10,
             "f" : 2.718,
             "s" : "Albuquerque"
             }
        plist = Teuchos.PyDictParameterList(d)
        for key in plist.dict():
            self.assertEqual(plist[key], plist.get(key))

    def test__len__(self):
        "Test Teuchos.PyDictParameterList __len__ method"
        d = {"i" : 10,
             "f" : 2.718,
             "s" : "Albuquerque",
             "d" : { "a" : 1, "b" : 2}
             }
        plist = Teuchos.PyDictParameterList(d)
        self.assertEqual(len(plist), len(d))

    def test__ne__1(self):
        "Test Teuchos.PyDictParameterList __ne__ method with dictionary"
        d1 = { "a" : 2 }
        d2 = { "a" : 2, "b" : 1 }
        plist = Teuchos.PyDictParameterList(d1)
        self.assertEqual(plist != d1, False)
        self.assertEqual(plist != d2, True )

    def test__ne__2(self):
        "Test Teuchos.PyDictParameterList __ne__ method with ParameterList"
        d = { "a" : 2, "b" : 1 }
        pdplist = Teuchos.PyDictParameterList(d)
        plist   = Teuchos.ParameterList(pdplist)
        self.assertEqual(pdplist != plist, False)
        plist.set("c",3)
        self.assertEqual(pdplist != plist, True )

    def test__repr__(self):
        "Test Teuchos.PyDictParameterList __repr__ method"
        plist = Teuchos.PyDictParameterList({"a":1,"b":2})
        d = eval(repr(plist))
        self.assertEqual(isinstance(d, dict), True)

    def test__setitem__1(self):
        "Test Teuchos.PyDictParameterList __setitem__ method"
        self.plist[self.name] = 2006
        self.assertEqual(self.plist.get(self.name), 2006)

    def test__setitem__2(self):
        "Test Teuchos.PyDictParameterList __setitem__ method for a ParameterList"
        plist = Teuchos.ParameterList()
        self.plist[self.name] = plist
        self.assertEqual(isinstance(self.plist.get(self.name),
                                    Teuchos.ParameterList), True)

    def test__setitem__3(self):
        "Test Teuchos.PyDictParameterList __setitem__ method for a PyDictParameterList"
        plist = Teuchos.PyDictParameterList()
        self.plist[self.name] = plist
        self.assertEqual(isinstance(self.plist.get(self.name),
                                    Teuchos.ParameterList), True)

    def test__str__(self):
        "Test Teuchos.PyDictParameterList __str__ method"
        plist = Teuchos.PyDictParameterList({"a":1,"b":2})
        d = eval(str(plist))
        self.assertEqual(isinstance(d, dict), True)

    def testDict(self):
        "Test Teuchos.PyDictParameterList dict method"
        d = {"maxits"         : 100,
             "tolerance"      : 1.0e-6,
             "preconditioner" : "diagonal"
             }
        plist = Teuchos.PyDictParameterList(d)
        self.assertEqual(plist.dict(), d)

    def testItems(self):
        "Test Teuchos.PyDictParameterList items method"
        d = {"maxits"         : 100,
             "tolerance"      : 1.0e-6,
             "preconditioner" : "diagonal"
             }
        plist = Teuchos.PyDictParameterList(d)
        self.assertEqual(plist.items(), d.items())

    def testKeys(self):
        "Test Teuchos.PyDictParameterList keys method"
        d = {"maxits"         : 100,
             "tolerance"      : 1.0e-6,
             "preconditioner" : "diagonal"
             }
        plist = Teuchos.PyDictParameterList(d)
        self.assertEqual(plist.keys(), d.keys())

    def testValues(self):
        "Test Teuchos.PyDictParameterList values method"
        d = {"maxits"         : 100,
             "tolerance"      : 1.0e-6,
             "preconditioner" : "diagonal"
             }
        plist = Teuchos.PyDictParameterList(d)
        self.assertEqual(plist.values(), d.values())

    def testUpdate1(self):
        "Test Teuchos.PyDictParameterList update method for a PyDictParameterList"
        intName    = "int parameter"
        intValue   = 8
        floatName  = "float parameter"
        floatValue = 3.14
        self.plist.set(intName,  intValue  )
        self.plist.set(floatName,floatValue)
        newList = Teuchos.PyDictParameterList()
        newList.update(self.plist)
        self.assertEqual(newList.get(intName  ), intValue  )
        self.assertEqual(newList.get(floatName), floatValue)
        self.assertEqual(newList.isSynchronized(), True)

    def testUpdate2(self):
        "Test Teuchos.PyDictParameterList update method for a ParameterList"
        plist = Teuchos.ParameterList()
        plist.set("int parameter"  , 8)
        plist.set("float parameter", 3.14)
        self.plist.update(plist)
        self.assertEqual(self.plist.isSynchronized(), True)

    def testUpdate3(self):
        "Test Teuchos.PyDictParameterList update method for a dictionary"
        d = {"int parameter"   : 8,
             "float parameter" : 3.14 }
        self.plist.update(d)
        self.assertEqual(self.plist.isSynchronized(), True)

    def test__iter__(self):
        "Test Teuchos.PyDictParameterList __iter__ method"
        values = range(10)
        for i in values:
            self.plist.set(str(i), i)
        keys = self.plist.keys()
        for p in self.plist:        # This invokes the __iter__ method
            self.assertEqual(p in keys, True)

####################################################################

if __name__ == "__main__":

    # Create the test suite object
    suite = unittest.TestSuite()

    # Add the test cases to the test suite
    suite.addTest(unittest.makeSuite(PyDictParameterListTestCase))

    # Create a communicator
    #comm    = Epetra.PyComm()
    #iAmRoot = comm.MyPID() == 0
    iAmRoot = True

    # Run the test suite
    if iAmRoot: print >>sys.stderr, \
       "\n***********************************\n" + \
       "Testing Teuchos.PyDictParameterList\n" + \
       "***********************************\n"
    verbosity = options.verbosity * int(iAmRoot)
    result = unittest.TextTestRunner(verbosity=verbosity).run(suite)

    # Exit with a code that indicates the total number of errors and failures
    #errsPlusFails = comm.SumAll(len(result.errors) + len(result.failures))
    errsPlusFails = len(result.errors) + len(result.failures)
    if errsPlusFails == 0 and iAmRoot: print "End Result: TEST PASSED"
    sys.exit(errsPlusFails)
