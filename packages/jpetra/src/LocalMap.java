/*
 * LocalMap.java
 *
 * Created on July 23, 2001, 2:14 PM
 */

package Jpetra;

/**
 *
 * @author  Michael William Boldt
 * @version 
 */
public class LocalMap extends Jpetra.Map {

    /** Creates new LocalMap */
    public LocalMap(int numProcessElements, int indexBase, Comm comm) {
	super(numProcessElements, numProcessElements, indexBase, comm);
	if(checkInput() != 0) {
	    System.out.println("Replicated Local Map not the same size onf all PEs");
	    System.exit(1);
	}
    }

    public LocalMap(LocalMap map) {
	super(map);
	isDistributedGlobal = false;
	if(checkInput() != 0) {
	    System.out.println("Replicated Local Map not the same size on all PEs");
	    System.exit(1);
	}
    }
    
    private int checkInput() {
	isDistributedGlobal = false;
	int [] tmp = new int [2];
	int [] res = new int [2];
	tmp[0] = numProcessElements;
	tmp[1] = - numProcessElements;
	res = comm.maxAll(tmp);

	int tmp1 = res[0];
	int tmp2 = - res[1];

	if(tmp1 == tmp2) return 0;
	else return -1;
    }
}
