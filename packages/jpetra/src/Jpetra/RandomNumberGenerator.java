package Jpetra;

/**
 *
 * @author  Jason Cross
 */
public interface RandomNumberGenerator {
    double getDouble();
    int getInt();
    void setRange(double min, double max);
}
