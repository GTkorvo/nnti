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
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************

# group classes by type/function in order to setup a correct build order
core =  Output.java \
	JpetraException.java \
	JpetraObject.java
misc_support =  Blas.java \
		NetlibBlas.java \
		JpetraBlas.java
ccj_support =   CcjSupport/CcjReduceAllSumIntArray.java  \
		CcjSupport/CcjReduceDoubleMaxArray.java \
		CcjSupport/CcjGatherDoubleArray.java \
		CcjSupport/CcjReduceDoubleMinArray.java \
		CcjSupport/CcjGatherIntArray.java \
		CcjSupport/CcjReduceIntMaxArray.java \
		CcjSupport/CcjGatherSerializableArray.java \
		CcjSupport/CcjReduceIntMinArray.java \
		CcjSupport/CcjScatterIntArray.java \
		CcjSupport/CcjScatter2DArray.java \
		CcjSupport/CcjReduceAllSumDoubleArray.java \
		CcjSupport/CcjReduceMaxSerializable.java \
		CcjSupport/CcjReduceMinSerializable.java \
		CcjSupport/CcjLink.java
comm =  Comm.java \
	Distributor.java \
	CcjComm.java \
	Directory.java \
	BasicDirectory.java \
	SerialComm.java  \
	DistObject.java \
        NullDistObject.java \
	Export.java \
	Import.java \
	CcjDistributor.java \
	SerialDistributor.java
util =  Time.java \
	Util.java \
        FlopCounter.java \
        RandomNumberGenerator.java \
        JavaRandom.java \
        RngPackRandom.java
data_structs =  ElementSpace.java \
		VectorSpace.java \
		Graph.java \
		CisMatrix.java \
		MultiVector.java \
		Vector.java
mmio =  MatrixMarketIO/VectorInfo.java \
	MatrixMarketIO/VectorSize.java \
	MatrixMarketIO/MatrixInfo.java \
	MatrixMarketIO/MatrixSize.java \
	MatrixMarketIO/MatrixVectorReader.java \
	MatrixMarketIO/CisMatrixReader.java \
	MatrixMarketIO/CisMatrixWriter.java \
	MatrixMarketIO/GraphReader.java \
	MatrixMarketIO/GraphWriter.java \
	MatrixMarketIO/MultiVectorReader.java \
	MatrixMarketIO/MultiVectorWriter.java 



#the compilation will fail if the classes are not compiled in the proper order		  
all:
	javac -d ../.. $(core) $(misc_support) $(util) $(ccj_support) $(comm) $(data_structs) $(mmio)
jar:
	jar -cf ../../Jpetra.jar ../../Jpetra/*.class ../../Jpetra/CcjSupport/*.class ../../Jpetra/MatrixMarketIO/*.class	
clean:
	rm -f ../../Jpetra/*.class ../../Jpetra/MatrixMarketIO/*.class ../../Jpetra/CcjSupport/*.class
