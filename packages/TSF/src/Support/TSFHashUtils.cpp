#include "TSFHashUtils.h"

using namespace TSF;


const int TSFHashUtils::primeCount_ = 29;
const int TSFHashUtils::primes_[] 
= {101, 163, 271, 443, 733, 1187, 1907, 3061, 
	 4919, 7759, 12379, 19543, 30841, 48487, 75989, 
	 119089, 185971, 290347, 452027, 703657, 1093237,
	 1695781, 2627993, 4067599, 6290467, 9718019, 
	 15000607, 23133937, 35650091};


int TSFHashUtils::nextPrime(int newCapacity) 
{
	if (newCapacity > primes_[primeCount_-1]) 
		TSFError::raise("TSFIntHashtable::nextPrime() overflow");

	for (int i=0; i<primeCount_; i++)
		{
			if (newCapacity <= primes_[i])
				{
					return primes_[i];
				}
		}
	// we shouldn't get here.
	TSFError::raise("TSFHashUtils::nextPrime() logic error");
	return 0;
}
