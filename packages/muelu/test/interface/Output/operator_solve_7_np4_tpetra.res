verbosity = high
coarse: max size = 100
transpose: use implicit = 1
max levels = 2
number of equations = 1   [default]
level 1 -> 
 P = Teuchos::RCP<Xpetra::Matrix<ignored> >{ptr=,node=,strong_count=6,weak_count=0}
 Nullspace = Teuchos::RCP<Xpetra::MultiVector<double, int, int, Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial> > >{ptr=,node=,strong_count=6,weak_count=0}

Clearing old data (if any)
MueLu::Amesos2Smoother: using "Superlu"
MueLu::AmesosSmoother: using "Superlu"
Using default factory (MueLu::SmootherFactory{pre = MueLu::DirectSolver{type = }, post = null}) for building 'CoarseSolver'.
Using default factory (MueLu::SmootherFactory{pre = MueLu::TrilinosSmoother{type = RELAXATION}, post = pre}) for building 'Smoother'.
Level 0
 Setup Smoother (MueLu::Ifpack2Smoother{type = RELAXATION})
  "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 1, damping factor: 1, Global matrix dimensions: [10000, 10000], Global nnz: 49600}
MueLu::Amesos2Smoother: using "Superlu"
MueLu::AmesosSmoother: using "Superlu"
Using default factory (MueLu::SmootherFactory{pre = MueLu::DirectSolver{type = }, post = null}) for building 'CoarseSolver'.
Using default factory (MueLu::SmootherFactory{pre = MueLu::TrilinosSmoother{type = RELAXATION}, post = pre}) for building 'Smoother'.
Level 1
 Computing Ac (MueLu::RAPFactory)
  MxM: A x P
   Matrix product nnz per row estimate = 5
  MxM: P' x (AP) (implicit)
   Matrix product nnz per row estimate = 5
  Ac size =  1700 x 1700, nnz = 15318
  Ac Load balancing info
  Ac   # active processes: 4/4
  Ac   # rows per proc   : avg = 4.25e+02,  dev =   0.0%,  min =   +0.0%,  max =   +0.0%
  Ac   #  nnz per proc   : avg = 3.83e+03,  dev =   0.0%,  min =   -0.0%,  max =   +0.0%
  Ac Communication info
  Ac   # num export send : avg =     0.00,  dev =   0.00,  min =    0.0 ,  max =    0.0
  Ac   # num import send : avg = 7.28e+01,  dev =   0.7%,  min =   -1.0%,  max =   +0.3%
  Ac   # num msgs        : avg =     3.00,  dev =   0.00,  min =    3.0 ,  max =    3.0
  Ac   # min msg size    : avg = 1.25e+00,  dev =  40.0%,  min =  -20.0%,  max =  +60.0%
  Ac   # max msg size    : avg = 3.75e+01,  dev =   1.5%,  min =   -1.3%,  max =   +1.3%
 Setup Smoother (MueLu::Amesos2Smoother{type = <ignored>})

--------------------------------------------------------------------------------
---                            Multigrid Summary                             ---
--------------------------------------------------------------------------------
Number of levels    = 2
Operator complexity = 1.31

matrix  rows    nnz  nnz/row procs
A 0    10000  49600     4.96  4
A 1     1700  15318     9.01  4

Smoother (level 0) both : "Ifpack2::Relaxation": {Initialized: true, Computed: true, Type: Symmetric Gauss-Seidel, sweeps: 1, damping factor: 1, Global matrix dimensions: [10000, 10000], Global nnz: 49600}

Smoother (level 1) pre  : <Direct> solver interface
Smoother (level 1) post : no smoother

