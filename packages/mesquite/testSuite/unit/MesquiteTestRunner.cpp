#include "MesquiteTestRunner.hpp"
#include "MsqTimer.hpp"
#include "cppunit/Test.h"
#include "cppunit/TestSuite.h"
#include "cppunit/TestResult.h"
#include "cppunit/TestFailure.h"
#include "cppunit/Exception.h"

#include <algorithm>

const unsigned char Mesquite::TestRunner::INDENT_SIZE = 2;
static bool last_test_succeeded = true;

void Mesquite::TestRunner::indent()
{
  unsigned int local_indent = indentLevel;
  while (local_indent--)
    *mOut << ' ';
}

void Mesquite::TestRunner::push_timer(Mesquite::Timer* timer)
{
  mTimers.push(timer);
}

Mesquite::Timer* Mesquite::TestRunner::pop_timer()
{
  Mesquite::Timer* rv = mTimers.top();
  mTimers.pop();
  return rv;
}

Mesquite::TestRunner::TestRunner() 
    : indentLevel(0)
{}

Mesquite::TestRunner::~TestRunner()
{
  delete_all_tests();
}

void Mesquite::TestRunner::delete_all_tests()
{
  for (std::vector<CppUnit::Test*>::iterator iter = mTests.begin();
       iter != mTests.end();
       ++iter)
  {
    delete *iter;
  }
  mTests.clear();
}

bool Mesquite::TestRunner::run(const std::string& name_of_run,
                               std::ostream& out_stream)
{
  mOut = &out_stream;
  *mOut << "Start of test run: " << name_of_run << "\n\n";
  
    // Run each test
  CppUnit::TestResult result;
  myResult = &result;
  result.addListener(this);
  for (std::vector<CppUnit::Test*>::iterator iter = mTests.begin();
       iter != mTests.end();
       ++iter)
    (*iter)->run(&result);
  
  *mOut << "\nEnd of test run: " << name_of_run << std::endl;

  return true;
}

void Mesquite::TestRunner::add_test(CppUnit::Test *test)
{
  mTests.push_back(test);
}

// This function is called by the TestRunner just before
// a test begins.
void Mesquite::TestRunner::startTest(CppUnit::Test *test)
{
  last_test_succeeded = true;
  
    // Indent
  indent();
  
    // Output a header
  *mOut << "Beginning of Test : " << test->getName();
  if (test->is_work_in_progress())
    *mOut << " (work in progress)";
  *mOut << std::endl;
  
    // increase the indent level
  indentLevel += TestRunner::INDENT_SIZE;

    // Start a timer
  push_timer(new Mesquite::Timer);
}

// This function is called by the TestRunner just before
// a suite begins.
void Mesquite::TestRunner::startSuite(CppUnit::TestSuite *suite)
{
    // Indent
  indent();
  
    // Output a header
  *mOut << "Beginning of Test Suite : " << suite->getName()
        << " (" << suite->countTestCases() << " tests)" << std::endl;

    // See if we've already run this suite
  if (suite->key() != 0 &&
      std::find(completedSuites.begin(),
                completedSuites.end(),
                suite->key()) != completedSuites.end())
  {
      // We've run this before...
    myResult->stop();
    indent();
    *mOut << "Tests run previously, skipping suite..." << std::endl;
  }
  else
  {
      // Add this test to the list of tests we've already started to run
    completedSuites.push_back(suite->key());
    
      // increase the indent level
    indentLevel += TestRunner::INDENT_SIZE;
    
      // Add a timer
    push_timer(new Mesquite::Timer);
      // Add a failure counter
    failureCounters.push(0);
  }
}

// This function is called if a test fails, either
// intentionally or unintentionally.
void Mesquite::TestRunner::addFailure(const CppUnit::TestFailure &failure)
{
  last_test_succeeded = false;
  if (!failure.failedTest()->is_work_in_progress())
    failureCounters.top() += 1;
  
    // Indicate whether error or failure.
    // An error is something you didn't specifically
    // look for or expect.  A failure is something you
    // looked for, but didn't like what you found.
  indent();
  if (failure.isError())
    *mOut << "***ERROR***\n";
  else
    *mOut << "*** FAILURE ***\n";
  
    // If it's an exception, say so.
  if (failure.thrownException())
  {
    indent();
    *mOut << "Error caught: " << failure.thrownException()->what()
          << std::endl;
  }

    // If we know where the error occurred, indicate it
  if (failure.sourceLine().isValid())
  {
    indent();
    *mOut << "Problem occured on line " << failure.sourceLine().lineNumber()
          << " of " << failure.sourceLine().fileName() << std::endl;
  }
}

// This function is called just after a test completes.
void Mesquite::TestRunner::endTest(CppUnit::Test *test)
{
    // Pop the timer
  Mesquite::Timer *timer = pop_timer();
  double elapsed_time = timer->since_birth();
  delete timer;
  
    // Decrease the indent level
  indentLevel -= TestRunner::INDENT_SIZE;
  
    // Output a footer
  indent();
  *mOut << test->getName();
  if (last_test_succeeded)
  {
    *mOut << " completed successfully in "
          << elapsed_time << " seconds" << std::endl;
  }
  else
  {
    *mOut << " failed after "<< elapsed_time
          << " seconds" << std::endl;
    if (test->is_work_in_progress())
    {
      indent();
      *mOut << "Test is work in progress, failure counted as success"
            << std::endl;
    }
  }
}

// This function is called just after a test completes.
void Mesquite::TestRunner::endSuite(CppUnit::TestSuite *suite)
{
  if (myResult->shouldStop())
  {
    myResult->reset();
    return;
  }
  
    // Pop the timer
  Mesquite::Timer *timer = pop_timer();
  double elapsed_time = timer->since_birth();
  delete timer;
  
    // Pop the success counter
  int failure_count = failureCounters.top();
  failureCounters.pop();
  if (!failureCounters.empty())
    failureCounters.top() += failure_count;
  
    // Decrease the indent level
  indentLevel -= TestRunner::INDENT_SIZE;
  
    // Output a footer
  indent();
  *mOut << suite->getName() << " Test Suite completed in "
        << elapsed_time << " seconds, ";
  if (failure_count)
  {
    *mOut << failure_count << " of " << suite->countTestCases()
          << " tests failed" << std::endl;
  }
  else
  {
    *mOut << "All " << suite->countTestCases() << " tests succeeded"
          << std::endl;
  }
}
