#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <vector>
#include <iostream>

#include <boost/timer.hpp>
#include <boost/foreach.hpp>

void force_calculation( long long sum )
{
  static int count = 0;
  //volatile to prevent the compiler from optimizing the loops away
  static volatile long long value;

  if( count && value != sum )
  {
    STKUNIT_EXPECT_TRUE(false);
  }

  value = sum;
  ++count;
}

struct sum_func
{
  sum_func() : sum(0) {}
  long long sum;
  void operator()( unsigned i )
  {
    sum += i;
  }
};

STKUNIT_UNIT_TEST( loop_iteration, loop_iteration)
{
  boost::timer  total_time;

  const size_t vector_size = 100000000;
  const size_t num_repetitions = 20;

  std::vector<unsigned> data;
  data.reserve(vector_size);

  {
    boost::timer timer;
    std::srand(10);
    std::generate_n(std::back_inserter(data), vector_size, &std::rand);
    double elapsedtime = timer.elapsed();
    std::cout << "Construction time: " << elapsedtime << " seconds." << std::endl;
  }
  {
    boost::timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      sum += std::accumulate( data.begin(), data.end(), 0ll );
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Iterator accumulate took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    boost::timer timer;
    unsigned* begin = &data.front();
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      sum += std::accumulate( begin, begin+data.size(), 0ll );
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Pointer accumulate took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    boost::timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      for( std::vector<unsigned>::iterator it = data.begin(), iend = data.end(); it != iend; ++it ) {
        sum += *it;
      }
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Iterator for loop took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    boost::timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      unsigned* it = &data.front();
      for( unsigned* iend = it+data.size(); it != iend; ++it ) {
        sum += *it;
      }
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Pointer for loop took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    boost::timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      for( std::size_t i = 0, size = data.size(); i != size; ++i ) {
        sum += data[i];
      }
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Index for (size_t) loop took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    boost::timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      for( int i = 0, size = data.size(); i != size; ++i ) {
        sum += data[i];
      }
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Index for (int) loop took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    boost::timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      sum += std::for_each( data.begin(), data.end(), sum_func() ).sum;
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Index for_each took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    boost::timer timer;
    unsigned* begin = &data.front();
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      sum += std::for_each( begin, begin+data.size(), sum_func() ).sum;
    }
    double elapsedtime = timer.elapsed();
    std::cout << "Pointer for_each took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }
  {
    boost::timer timer;
    long long sum = 0;
    for (size_t n = 0; n < num_repetitions; ++n) {
      BOOST_FOREACH( unsigned i, data ) {
        sum += i;
      }
    }
    double elapsedtime = timer.elapsed();
    std::cout << "BOOST_FOREACH took " << elapsedtime << " seconds." << std::endl;
    force_calculation(sum);
  }

  double elapsedtime = total_time.elapsed();

  std::cout << "Total time: " << elapsedtime << " seconds." << std::endl;
}
