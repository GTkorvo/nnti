//
//  File : BlockArnoldiPetraExGen.cpp
//
//  This example computes the eigenvalues of smallest magnitude of the discretized 1D Laplacian
//  equation using the block Implicitly-Restarted Arnoldi method.  This problem shows the construction
//  of an inner-outer iteration using Belos as the linear solver within Anasazi.  An Ifpack preconditioner
//  is constructed to precondition the linear solver.  This operator is discretized using finite elements
//  and constructed as an Epetra matrix, then passed into the AnasaziPetraMat to be used in the construction
//  of the Krylov decomposition.  The specifics of the block Arnoldi method can be set by the user.

#include "AnasaziPetraInterface.hpp"
#include "BelosPetraInterface.hpp"
#include "AnasaziBlockArnoldi.hpp"
#include "AnasaziCommon.hpp"
#include "Ifpack_CrsIct.h"
#include "Epetra_CrsMatrix.h"
#include "BelosEpetraOperator.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[]) {
	int ierr = 0, i, j;

#ifdef EPETRA_MPI

	// Initialize MPI

	MPI_Init(&argc,&argv);

	int size, rank; // Number of MPI processes, My process ID

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

	int size = 1; // Serial case (not using MPI)
	int rank = 0;

#endif


#ifdef EPETRA_MPI
	Epetra_MpiComm & Comm = *new Epetra_MpiComm(MPI_COMM_WORLD);
#else
	Epetra_SerialComm & Comm = *new Epetra_SerialComm();
#endif

	int MyPID = Comm.MyPID();
	int NumProc = Comm.NumProc();
	cout << "Processor "<<MyPID<<" of "<< NumProc << " is alive."<<endl;

	bool verbose = (MyPID==0);

        //  Dimension of the matrix
        int NumGlobalElements = 1000;

        // We will use zero based indices
        int IndexBase = 0;

        // Construct a Map that puts approximately the same number of
        // equations on each processor.

        Epetra_Map& Map = *new Epetra_Map(NumGlobalElements, 0, Comm);

        // Get update list and number of local equations from newly created Map.
        
        int NumMyElements = Map.NumMyElements();

        int * MyGlobalElements = new int[NumMyElements];
        Map.MyGlobalElements(MyGlobalElements);

        // Create an integer vector NumNz that is used to build the Petra Matrix.
        // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
        // on this processor
        int * NumNz = new int[NumMyElements];

        // We are building two tridiagonal matrices
        // So we need 2 off-diagonal terms (except for the first and last equation)

        for (i=0; i<NumMyElements; i++) {
                if (MyGlobalElements[i]==0 || MyGlobalElements[i] == NumGlobalElements-1)
                        NumNz[i] = 2;
                else
                        NumNz[i] = 3;  
        }

        // Create both the stiffness and mass Epetra_Matrix        
        Epetra_CrsMatrix& A = *new Epetra_CrsMatrix(Copy, Map, NumNz);
        Epetra_CrsMatrix& B = *new Epetra_CrsMatrix(Copy, Map, NumNz);

        const double one = 1.0;
        double *ValuesA = new double[2];
        double *ValuesB = new double[2];

	// Set values of stiffness matrix.
        double h = one /(NumGlobalElements + one);
        ValuesA[0] = -one/h; ValuesA[1] = -one/h;
        double diagA = 2.0/h;

	// Set values of mass matrix.
        h = one /(6.0*(NumGlobalElements + one));
        ValuesB[0] = one/h; ValuesB[1] = one/h;
        double diagB = 4.0/h;
        int *Indices = new int[2];
        int NumEntries;

        for (i=0; i<NumMyElements; i++)
        {
                if (MyGlobalElements[i]==0)
                {
                        Indices[0] = 1;
                        NumEntries = 1;
                        assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, ValuesA+1, Indices)==0);
                        assert(B.InsertGlobalValues(MyGlobalElements[i], NumEntries, ValuesB+1, Indices)==0);
                }
                else if (MyGlobalElements[i] == NumGlobalElements-1)
                {
                        Indices[0] = NumGlobalElements-2;
                        NumEntries = 1;
                        assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, ValuesA, Indices)==0);
                        assert(B.InsertGlobalValues(MyGlobalElements[i], NumEntries, ValuesB, Indices)==0);
                }
                else
                {
                        Indices[0] = MyGlobalElements[i]-1;
                        Indices[1] = MyGlobalElements[i]+1;
                        NumEntries = 2;
                        assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, ValuesA, Indices)==0);
                        assert(B.InsertGlobalValues(MyGlobalElements[i], NumEntries, ValuesB, Indices)==0);
                }
                // Put in the diagonal entry
                assert(A.InsertGlobalValues(MyGlobalElements[i], 1, &diagA, MyGlobalElements+i)==0);
                assert(B.InsertGlobalValues(MyGlobalElements[i], 1, &diagB, MyGlobalElements+i)==0);
        }
         
        // Finish up
        assert(A.TransformToLocal()==0);
        A.SetTracebackMode(1); // Shutdown Epetra Warning tracebacks
        assert(B.TransformToLocal()==0);
        B.SetTracebackMode(1); // Shutdown Epetra Warning tracebacks


	//
        //*****Select the Preconditioner*****
        //
        if (verbose) cout << endl << endl;
        if (verbose) cout << "Constructing ICT preconditioner" << endl;
        int Lfill = 0;
        if (argc > 2) Lfill = atoi(argv[2]);
        if (verbose) cout << "Using Lfill = " << Lfill << endl;
        int Overlap = 0;
        if (argc > 3) Overlap = atoi(argv[3]);
        if (verbose) cout << "Using Level Overlap = " << Overlap << endl;
        double Athresh = 0.0;
        if (argc > 4) Athresh = atof(argv[4]);
        if (verbose) cout << "Using Absolute Threshold Value of " << Athresh << endl;
        double Rthresh = 1.0;
        if (argc >5) Rthresh = atof(argv[5]);
        if (verbose) cout << "Using Relative Threshold Value of " << Rthresh << endl;
        double dropTol = 1.0e-6;
        //
        Ifpack_CrsIct* ICT = 0;
        //
        if (Lfill > -1) {
                ICT = new Ifpack_CrsIct(A, dropTol, Lfill);

                ICT->SetAbsoluteThreshold(Athresh);
                ICT->SetRelativeThreshold(Rthresh);
                int initerr = ICT->InitValues(A);
                if (initerr != 0) cout << "InitValues error = " << initerr;
                assert(ICT->Factor() == 0);
        } 
        //
        bool transA = false;
        double Cond_Est;
        ICT->Condest(transA, Cond_Est);
        if (verbose) {
                cout << "Condition number estimate for this preconditoner = " << Cond_Est << endl;
                cout << endl;
        } 
        Epetra_Operator& prec = dynamic_cast<Epetra_Operator&>(*ICT);
        //
        // call the ctor for the preconditioning object
        //
	BelosPetraMat<double> BelosMat(A);
        BelosPetraPrec<double> EpetraOpPrec(prec);

	//*******************************************************
	// Set up Belos Block GMRES operator for inner iteration
	//*******************************************************

        int block = 3;  // blocksize used by solver
	int numrhs = block; // number of right hand sides
        int maxits = NumGlobalElements - 1; // maximum number of iterations to run
        double btol = 1.0e-10;  // relative residual tolerance

	BelosEpetraOperator<double> BelosOp( BelosMat, EpetraOpPrec, "BlockCG", btol, maxits, block, 0, 0 );

	//************************************
	// Start the block Arnoldi iteration
	//************************************
	//
	//  Variables used for the Block Arnoldi Method
	//
	int length = 20;
	int nev = 10;
	double tol = 1.0e-8;
	string which="LM";
	int step = 5;
	int restarts = 5;

	// create a PetraAnasaziVec. Note that the decision to make a view or
	// or copy is determined by the petra constructor called by AnasaziPetraVec.
	// This is possible because I pass in arguements needed by petra.

	AnasaziPetraVec<double> ivec(Map, block);
	ivec.MvRandom();

	// call the ctor that calls the petra ctor for a matrix

	AnasaziPetraMat<double> Amat(A);
	AnasaziPetraMat<double> Bmat(B);
	AnasaziPetraGenOp<double> Aop(BelosOp, B);	
	AnasaziEigenproblem<double> MyProblem(&Amat, &Bmat, &Aop, &ivec);

	// initialize the Block Arnoldi solver
	Anasazi::BlockArnoldi<double> MyBlockArnoldi(MyProblem, tol, nev, length, block, 
						which, step, restarts);
	
	// inform the solver that the problem is symmetric
	//MyBlockArnoldi.setSymmetric(true);
	MyBlockArnoldi.setDebugLevel(0);

#ifdef UNIX
	Epetra_Time & timer = *new Epetra_Time(Comm);
#endif

	// iterate a few steps (if you wish)
	//MyBlockArnoldi.iterate(5);

	// solve the problem to the specified tolerances or length
	MyBlockArnoldi.solve();

#ifdef UNIX
	double elapsed_time = timer.ElapsedTime();
	double total_flops = A.Flops();
	double MFLOPs = total_flops/elapsed_time/1000000.0;
#endif

	// obtain results directly
	double* resids = MyBlockArnoldi.getResiduals();
	double* evalr = MyBlockArnoldi.getEvals(); 

	// retrieve real and imaginary parts of the eigenvectors
	AnasaziPetraVec<double> evecr(Map, nev);
	MyBlockArnoldi.getEvecs( evecr );

	AnasaziDenseMatrix<double> dmatr(nev,nev);
	MyProblem.AInProd( one, evecr, evecr, dmatr );
	double* ptr_dmatr = dmatr.getarray();
	double compeval;

	cout<<"Actual Eigenvalues (obtained by Rayleigh quotient) : "<<endl;
	cout<<"Real Part \t Rayleigh Error"<<endl;
	for (i=0; i<nev; i++) {
		compeval = ptr_dmatr[i*nev+i];
		cout<<compeval<<"\t"<<abs(compeval-one/evalr[i])<<endl;
	}

	// output results to screen
	MyBlockArnoldi.currentStatus();


#ifdef UNIX
	if (verbose)
		cout << "\n\nTotal MFLOPs for Arnoldi = " << MFLOPs << " Elapsed Time = "<<  elapsed_time << endl;
#endif


	// Release all objects
	delete [] resids, evalr;

	delete &A, &B;
	delete &Map;
	delete &Comm;

	return 0;
}
