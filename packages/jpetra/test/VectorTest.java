/*
 * VectorTest.java
 *
 * Created on February 9, 2004, 4:27 PM
 */

package test;

import Jpetra.*;
import java.util.Random;

/**
 *
 * @author  jc
 */
public class VectorTest extends JpetraObject {
    
    /** Creates a new instance of VectorTest */
    public VectorTest(String[] args) {
        this.initializeOutput();
        if(args.length > 0 && args[0].equals("-v")) setRootPrint("VERBOSE", true);
        if(args.length > 0 && args[0].equals("-d")) {
            setRootPrint("DEBUG", true);
            setRootPrint("VERBOSE", true);
        }
        
        int ierr = 0;
	println("STD", "Starting VectorTest...");
	ierr += unitTests();

	// finish up
	if(ierr == 0)
		println("STD", "Vector test successful.");
	else
		println("STD", "Vector test failed.  Error #: " + ierr);
	
        
        
    }
    
    public int unitTests() {
	int ierr = 0;
	int returnierr = 0;
        SerialPlatform platformE = new SerialPlatform();
	SerialPlatform platformV = new SerialPlatform();

	//if(verbose) cout << "Starting unit tests for Vector<" << OTName << "," << STName << ">." << endl;

	//
	// code coverage section - just call functions, no testing
	//
	println("VERBOSE", "Starting code coverage section...");
	// constructors
	println("VERBOSE", "Constructors...");
	// taking a VectorSpace
	int ESlength = 10;
	int ESindexBase = 2;
	ElementSpace elementspace = new ElementSpace(ESlength, ESindexBase, platformE);
	VectorSpace vectorspace = new VectorSpace(elementspace, platformV);
	Vector vector = new Vector(vectorspace);
	// taking a VectorSpace and a user array of entries
	double[] scalarArray = new double[ESlength];
        
        Random random = new Random();
	for(int i = 0; i < ESlength; i++)
		scalarArray[i] = random.nextDouble();
	Vector vector1a = new Vector(scalarArray, ESlength, vectorspace);
	// cpy ctr
	Vector v2 = new Vector(vector);

	// print
		//cout << "Overloaded << operator..." << endl;
		//cout << vector << endl;
            println("VERBOSE", "print method test:");
            vector.print(System.out);

	// attribute access
	println("VERBOSE", "Attribute access methods...");
	int temp = 0;
	temp = vector.getNumGlobalEntries();
	temp = vector.getNumMyEntries();

	// element access
	println("VERBOSE", "Element access methods...");
	double temp1 = vector.getIndex(1);
	vector.setIndex(0, temp1);

	println("VERBOSE", "Code coverage section finished.");

	//
	// actual testing section - affects return code
	//

	println("VERBOSE", "Starting actual testing section...");

	// default ctr initializing to zero
	println("VERBOSE", "Checking to see that default constructor initializes to zeros... ");
	Vector testVector1 = new Vector (vectorspace);
	int length = testVector1.getNumMyEntries();
	for(int i = 0; i < length; i++)
		if(testVector1.getIndex(i) != 0) {
			println("DEBUG", "element " + i + " = " + testVector1.getIndex(i) + ", should be zero");
			ierr++;
		}
		if(ierr == 0) 
			println("VERBOSE", "Passed");
		else
			println("VERBOSE", "Failed");
	returnierr += ierr;
	ierr = 0;

	// user array ctr initializing correctly
	println("VERBOSE", "Checking to see that user array constructor initializes values correctly... ");
	length = vector1a.getNumMyEntries();
	for(int i = 0; i < length; i++)
		if(vector1a.getIndex(i) != scalarArray[i])
                    
		if(ierr == 0) 
			println("VERBOSE","Passed");
		else
			println("VERBOSE","Failed");
	returnierr += ierr;
	ierr = 0;
	
	// changing data
	println("VERBOSE", "Changing data... ");
	double value2 = random.nextDouble();
	double value5 = random.nextDouble();
	double value0 = random.nextDouble();
	testVector1.setIndex(2, value2);
	testVector1.setIndex(5, value5);
	testVector1.setIndex(0, value0);
        
	if(testVector1.getIndex(2) != value2) {
		println("DEBUG", "element 2 = " + testVector1.getIndex(2) + ", should be " + value2);
		ierr++;
	}
	if(testVector1.getIndex(5) != value5) {
		println("DEBUG", "element 5 = " + testVector1.getIndex(5) + ", should be " + value5);
		ierr++;
	}
	if(testVector1.getIndex(0) != value0) {
		println("DEBUG", "element 0 = " + testVector1.getIndex(0) + ", should be " + value0);
		ierr++;
	}
		if(ierr == 0) 
			println("VERBOSE", "Passed");
		else
			println("VERBOSE", "Failed");
	returnierr += ierr;
	ierr = 0;

	// finish up
		if(returnierr == 0)
			println("VERBOSE", "VectorTest passed.");
		else
			println("VERBOSE","VectorTest failed.");
	return(returnierr);
    }
    
    public static void main(String[] args) {
        new VectorTest(args);
    }
}
