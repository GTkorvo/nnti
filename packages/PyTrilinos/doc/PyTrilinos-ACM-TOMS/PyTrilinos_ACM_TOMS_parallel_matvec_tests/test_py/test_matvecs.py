#! /usr/bin/env python
from PyTrilinos import Epetra, Triutils

Epetra.Init()
Comm = Epetra.PyComm()
nx = 4000 
ny = 1000 * Comm.NumProc()
Gallery = Triutils.CrsMatrixGallery("laplace_2d", Comm)
Gallery.Set("nx", nx)
Gallery.Set("ny", ny)
Gallery.Set("problem_size", nx * ny)
Gallery.Set("map_type", "linear")
Matrix = Gallery.GetMatrix()
LHS = Gallery.GetStartingSolution()
RHS = Gallery.GetRHS()
Time = Epetra.Time(Comm)

for i in xrange(10):
  Matrix.Multiply(False, LHS, RHS)

print Time.ElapsedTime()
Epetra.Finalize()
