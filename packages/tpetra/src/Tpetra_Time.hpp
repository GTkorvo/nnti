#ifndef _TPETRA_TIME_HPP_
#define _TPETRA_TIME_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_Comm.hpp"
#include <sys/time.h>
#include <sys/resource.h>

namespace Tpetra {

//! Tpetra::Time class, utility for timing code segments.
/*! To time a section of code, place it in between calls to resetStartTime and elapsedTime. 
   Note that resetStartTime is called by the constructor, so it is also possible to time a section of code 
   by placing it in between the creation of the Time object, and a call to elapsedTime. 
   A call to resetStartTime must be used for all subsequent timings though.

   A Comm object is also required to use Time, although I don't know why.
*/

class Time : public Tpetra::Object {

public:
  //! Default constructor
  Time(const Comm<double, int>& Comm);
  
  //! Copy Constructor
  Time(const Time& Time);
  
  //! Destructor
  virtual ~Time() {};
  
  //! returns current wall-clock time in seconds
  double wallTime() const;
  
  //! resets the timer to the current walltime
  void resetStartTime();
  
  //! returns the elapsed time in between when the timer was set, and the current walltime.
  double elapsedTime() const;


private:
  double startTime_;
  const Comm<double, int> * Comm_;
};

} // namespace Tpetra

#include "Tpetra_Time.cpp"

#endif // TPETRA_TIME_HPP_
