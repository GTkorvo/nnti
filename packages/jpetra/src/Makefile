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

# classes that currently don't compile or are not used
# Histogram.java VecView.java Graph.java

# group classes by type/function in order to setup a correct build order
core_classes = Output.java JpetraException.java JpetraObject.java
ccj_classes =  CcjComm.java CcjReduceAllSumIntArray.java CcjDistributor.java \
		CcjReduceDoubleMaxArray.java CcjGatherDoubleArray.java \
		CcjReduceDoubleMinArray.java CcjGatherIntArray.java \
		CcjReduceIntMaxArray.java CcjGatherSerializableArray.java \
		CcjReduceIntMinArray.java CcjLink.java CcjScatter2DArray.java \
		CcjReduceAllSumDoubleArray.java
comm_classes = Comm.java SerialComm.java Directory.java Blockmap.java Map.java \
		LocalMap.java Exporter.java Importer.java \
		DistObject.java Distributor.java SerialDistributor.java $(ccj_classes)
util_classes = Time.java Util.java
data_structs = ElementSpace.java CrsMatrix.java SerialDenseMatrix.java IntSerialDenseMatrix.java \
		SerialSymDenseMatrix.java SerialDenseVector.java IntSerialDenseVector.java \
		VectorSpace.java Vector.java MultiVector.java
		  
all:
	#the compilation will fail if the classes are not compiled in the proper order
	javac -d .. $(core_classes) $(util_classes) $(comm_classes) $(data_structs)
	
	# old build command
	#javac -d ..  BlockMap.java Comm.java CrsMatrix.java Directory.java Exporter.java Graph.java \
        #          Importer.java JpetraException.java JpetraObject.java LocalMap.java \
        #          Map.java MatView.java MultiVector.java SerialComm.java Time.java Util.java \
        #          Vector.java CcjGatherIntArray.java CcjGatherDoubleArray.java CcjReduceAllSumDoubleArray.java \
        #          CcjReduceAllSumIntArray.java CcjLink.java


clean:
	rm *~
