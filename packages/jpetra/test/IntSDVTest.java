import Jpetra.*;

public class IntSDVTest {
    public static void main(String[] args) {
        new IntSDVTest();
    }

    public IntSDVTest() {
        int length = 5;
        IntSerialDenseVector isdv = new IntSerialDenseVector(length);
        
        // fill up isdv
        int i;
        int k = 0;
        for(i=0; i < length; i++) {
            isdv.setElement(i, k++);
        }
        
        System.out.println("isdv:");
        isdv.print();
        
        IntSerialDenseVector isdv1 = new IntSerialDenseVector(isdv);

        System.out.println("isdv1:");        
        isdv1.print();
        
        System.out.println("isdm OneNorm: " + isdv.getOneNorm());
        System.out.println("isdv1 OneNorm: " + isdv1.getOneNorm());
        
        System.out.println("isdv InfoNorm: " + isdv.getInfNorm());
        System.out.println("isdv1 InfoNorm: " + isdv1.getInfNorm());
        
        isdv1.resize(length+2);

        System.out.println("isdv:");        
        isdv.print();
        System.out.println("isdv1:");        
        isdv1.print();
        
        System.out.println("isdv OneNorm: " + isdv.getOneNorm());
        System.out.println("isdv1 OneNorm: " + isdv1.getOneNorm());
        
        System.out.println("isdv InfoNorm: " + isdv.getInfNorm());
        System.out.println("isdv1 InfoNorm: " + isdv1.getInfNorm());
        
        isdv1.resize(length-2);
        isdv.size(length-3);
        
        System.out.println("isdv:");        
        isdv.print();
        System.out.println("isdv1:");        
        isdv1.print();
        
        System.out.println("isdv OneNorm: " + isdv.getOneNorm());
        System.out.println("isdv1 OneNorm: " + isdv1.getOneNorm());
        
        System.out.println("isdv InfoNorm: " + isdv.getInfNorm());
        System.out.println("isdv1 InfoNorm: " + isdv1.getInfNorm());
    }
}