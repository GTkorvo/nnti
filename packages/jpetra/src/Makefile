# ************************************************************************
# 
#               Java Implementation of the Petra Library
#                 Copyright (2004) Sandia Corporation
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

# group classes by type/function in order to setup a correct build order
core = Output.java JpetraException.java JpetraObject.java
misc_support = Blas.java NetlibBlas.java JpetraTreeMap.java
ccj_support = CcjReduceAllSumIntArray.java  \
		CcjReduceDoubleMaxArray.java CcjGatherDoubleArray.java \
		CcjReduceDoubleMinArray.java CcjGatherIntArray.java \
		CcjReduceIntMaxArray.java CcjGatherSerializableArray.java \
		CcjReduceIntMinArray.java CcjLink.java CcjScatter2DArray.java \
		CcjReduceAllSumDoubleArray.java
comm = Comm.java Distributor.java CcjComm.java Directory.java SerialComm.java  \
		 DistObject.java Export.java Import.java \
		 CcjDistributor.java SerialDistributor.java
util = Time.java Util.java
data_structs = ElementSpace.java VectorSpace.java Graph.java CisMatrix.java MultiVector.java
		  
all:
	#the compilation will fail if the classes are not compiled in the proper order
	javac -d .. $(core) $(misc_support) $(util) $(ccj_support) $(comm) $(data_structs)
	
	# old build command
	#javac -d ..  BlockMap.java Comm.java CrsMatrix.java Directory.java Exporter.java Graph.java \
        #          Importer.java JpetraException.java JpetraObject.java LocalMap.java \
        #          Map.java MatView.java MultiVector.java SerialComm.java Time.java Util.java \
        #          Vector.java CcjGatherIntArray.java CcjGatherDoubleArray.java CcjReduceAllSumDoubleArray.java \
        #          CcjReduceAllSumIntArray.java CcjLink.java


clean:
	rm *~
