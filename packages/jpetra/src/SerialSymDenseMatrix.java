package Jpetra;

public class SerialSymDenseMatrix extends SerialDenseMatrix {
    boolean upper;
    String uplo;
    
    public SerialSymDenseMatrix() {
        super();

        upper = false;
        uplo = "L";
    }
    
    public SerialSymDenseMatrix(double[][] matrixA, int lda, int numRowsCols) {
        super(matrixA, lda, numRowsCols);
        
        upper = false;
        uplo = "L";
    }
    
    public SerialSymDenseMatrix(SerialSymDenseMatrix source) {
        super(source);

        this.upper = source.getUpper();
        this.uplo = source.getUpLo();
    }
    
    public void CopyUpLoMatrix(boolean upper, double[] myA, int lda, int numRows) {   
        int i, j;        
        int offSetA, offSetB;
        
        if (upper) {
            for (j=1; j<numRows; j++) {
                offSetA = j;
                offSetB = j*lda;
                for (i=0; i<j; i++) {
                    myA[offSetA] = myA[offSetB++];
                    offSetA += lda;
                }
            }
        }
        else {
            for (i=1; i<numRows; i++) {
                offSetA = i;
                offSetB = i*lda;
                for (j=0; j<i; j++) {
                    myA[offSetB++] = myA[offSetA];
                    offSetA += lda;
                }
            }
        }
    }
    
    public double getInfNorm() {
        
        int i, j;
        
        double anorm = 0.0;
        int offSet;
        
        if (!this.upper) {
            for (j=0; j< this.numCols; j++) {
                double sum = 0.0;
                offSet = j + j*this.lda;
                for (i=j; i< this.numCols; i++) {
                    sum += Math.abs(this.matrixA[offSet++]);
                }
                offSet = j;
                for (i=0; i<j; i++) {
                    sum += Math.abs(this.matrixA[offSet]);
                    offSet += this.lda;
              }
              anorm = Math.max(anorm, sum);
            }
        }
        else {
            for (j=0; j< this.numCols; j++) {
                double sum = 0.0;
                offSet = j*this.lda;
                for (i=0; i<j; i++) {
                    sum += Math.abs(this.matrixA[offSet++]);
                }
                offSet = j + j*this.lda;
                for (i=j; i< this.numCols; i++) {
                    sum += Math.abs(this.matrixA[offSet]);
                    offSet += this.lda;
                }
              anorm = Math.max(anorm, sum);
            }
        }
        
        /*UpdateFlops(N_*N_);*/
        
        return(anorm);
    }
    
    public int scale(double scalarA) {
        int i, j;
        
        int offSet;
        
        if (!this.upper) {
            for (j=0; j< this.numCols; j++) {
                offSet = j + j*this.lda;
                for (i=j; i< this.numCols; i++) {
                    this.matrixA[offSet++] *= scalarA;
                }
            }
        }
        else {
            for (j=0; j< this.numCols; j++) {
                offSet =  j*this.lda;
                for (i=0; i<j; i++) {
                    matrixA[offSet++] *= scalarA;
                }
            }
        }
        
        /*UpdateFlops(N_*(N_+1)/2);*/
        
        return(0);
    }

    
    public double getOneNorm() {
        return(getInfNorm());
    }
    
    public boolean getUpper() {
        return this.upper;
    }
    
    public String getUpLo() {
        return uplo;
    }
}