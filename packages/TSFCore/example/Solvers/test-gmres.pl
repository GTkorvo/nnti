#!/usr/bin/perl -w
use strict;
use strict 'refs';
my $success = 1;
my $result;
$result = system ('./TSFCoreSolversGmresTest.exe ../../../../../packages/belos/test/BlockGmres/orsirr1.hb 500 1e-10 -v');
$success = 0 if (!$result);
$result = system ('./TSFCoreSolversGmresTest.exe ../../../../../packages/belos/test/BlockGmres/fidap036.hb 1500 1e-10 -v');
$success = 0 if (!$result);
$result = system ('./TSFCoreSolversGmresTest.exe ../../../../../packages/belos/test/BlockCG/bcsstk14.hb 1200 1e-10 -v');
$success = 0 if (!$result);
$result = system ('./TSFCoreSolversGmresTest.exe ../../../../../packages/belos/test/BlockCG/bcsstk15.hb 2100 1e-10 -v');
$success = 0 if (!$result);
if($success) {
  printf "\nFinal: Congraduations! all of the linear systems where solved the the given tolerances!\n";
}
else {
  printf "\nFinal: Congraduations! all of the linear systems where solved the the given tolerances!\n";
}
exit ($success ? 0 : -1 );
