function java = copyMatlabMatrix(mat)
import java.lang.Double;
mySize = size(mat);
m = mySize(1);
n = mySize(2);

java(m,n) = Double(0);

for tmpM = 1:m
    for tmpN = 1:n
        java(tmpM,tmpN) = Double(mat(tmpM,tmpN));
    end
end
