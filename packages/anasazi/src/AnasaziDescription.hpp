/*!
   \mainpage %Anasazi: A Block Eigensolvers Package

   \section intro Introduction

   %Anasazi is an extensible and interoperable framework for large-scale eigenvalue algorithms.
The motivation for this framework is to provide a generic interface to a collection of algorithms for 
solving large-scale eigenvalue problems.
Anasazi is interoperable because both the matrix and vectors (defining the
eigenspace) are considered to be opaque objects---only knowledge of the matrix and
vectors via elementary operations is necessary. An implementation of Anasazi
is accomplished via the use of interfaces. Current interfaces available include
Epetra and so any libraries that understand Epetra matrices and vectors (such
as AztecOO) may also be used in conjunction with Anasazi.

One of the goals of Anasazi is to allow the user the flexibility to specify the
data representation for the matrix and vectors and so leverage any existing software
investment. The algorithm that is currently available through Anasazi is block implicitly restarted Arnoldi. 

   \section contributors Anasazi Contributors

   The following people have contributed to the development of %Anasazi:

   <ul>
 	<li> Rich Lehoucq, Sandia National Labs, rblehou@sandia.gov
	<li> Heidi Thornquist, Sandia National Labs, hkthorn@sandia.gov
   </ul>

*/
