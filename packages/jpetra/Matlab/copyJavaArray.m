function mat = copyJavaArray(java)
mySize = size(java);
m = mySize(1);
mySize = size(java(1));
n = mySize(1);

mat = zeros(m,n);

for tmpM = 1:m
    for tmpN = 1:n
        mat(tmpM,tmpN) = double(java(tmpM,tmpN));
    end
end
