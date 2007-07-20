#!/bin/csh

echo "Running ARPACK examples"
cd ARPACK
rm ar.txt
foreach nx (50 100 250 500 1000)
   foreach ncv (50 100 150)
      ./nrun.csh $nx $ncv >> ar.txt
   end
end
cd ..

echo "Running Anasazi examples"
cd Anasazi
rm an.txt
foreach nx (50 100 250 500 1000)
   foreach ncv (49 99 149)
      ./nrun.csh $nx $ncv >> an.txt
   end
end
cd ..

echo "Done."
