package Jpetra;

import java.util.Random;

/**
 *
 * @author  Jason Cross
 */
public class JavaRandom extends JpetraObject implements RandomNumberGenerator {
    Random random;
    double min;
    double max;
    double range;
    
    /** Creates a new instance of JavaRandom */
    public JavaRandom() {
        this.random = new Random();
    }
    
    public JavaRandom(long seed) {
        this.random = new Random(seed);
    }
    
    public void setRange(double min, double max) {
        this.min = min;
        this.max = max;
        this.range = max - min;
    }
    
    public double getDouble() {
        return (random.nextDouble() * range) - min;
    }
    
    public int getInt() {
        return (int) random.nextInt((int) range) + (int) min;
    }
    
    public Random getRandom() {
        return this.random;
    }
}
