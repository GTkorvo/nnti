package Jpetra;

import edu.cornell.lassp.houle.RngPack.*;

/**
 *
 * @author  Jason Cross
 */
public class RngPackRandom extends JpetraObject implements RandomNumberGenerator {
    public static final int JAVA = 0;
    public static final int MT = 1;
    public static final int ECU = 2;
    public static final int LUX = 3;
    public static final int MAR = 3;
    
    double min;
    double max;    
    RandomElement random;

    public RngPackRandom(int generatorType, double min, double max) {
        this.min = min;
        this.max = max;
        
        if (generatorType == JAVA) {
            this.random = new RandomJava();
        }
        else if (generatorType == MT) {
            this.random = new RanMT();
        }
        else if (generatorType == ECU) {
            this.random = new Ranecu();
        }
        else if (generatorType == LUX) {
            this.random = new Ranlux();
        }
        else if (generatorType == MAR) {
            this.random = new Ranmar();
        }
    }
    
    public double getDouble() {
        return random.uniform(min, max);
    }
    
    public int getInt() {
        return random.choose((int) min, (int) max);
    }
    
    public void setRange(double min, double max) {
        this.min = min;
        this.max = max;
    }
    
    public RandomElement getRandomElement() {
        return this.random;
    }
    
}
