//  Copyright Andrey Semashev 2023.
//
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

//  This file is a copy of boost/timer.hpp from Boost.Timer v1, which
//  was deprecated and slated for removal. We are not using Boost.Timer v2
//  components to avoid having to link with its binary.
//
//  See http://www.boost.org/libs/timer for documentation.

#ifndef BOOST_PARAMETER_TEST_TIMER_HPP
#define BOOST_PARAMETER_TEST_TIMER_HPP

#include <ctime>
#include <limits>

namespace test {

//  timer  -------------------------------------------------------------------//

//  A timer object measures elapsed time.

//  It is recommended that implementations measure wall clock rather than CPU
//  time since the intended use is performance measurement on systems where
//  total elapsed time is more important than just process or CPU time.

//  Warnings: The maximum measurable elapsed time may well be only 596.5+ hours
//  due to implementation limitations.  The accuracy of timings depends on the
//  accuracy of timing information provided by the underlying platform, and
//  this varies a great deal from platform to platform.

class timer
{
 public:
  timer() { _start_time = std::clock(); } // postcondition: elapsed()==0
  void restart() { _start_time = std::clock(); } // post: elapsed()==0
  double elapsed() const                  // return elapsed time in seconds
    { return  double(std::clock() - _start_time) / CLOCKS_PER_SEC; }

  double elapsed_max() const   // return estimated maximum value for elapsed()
  // Portability warning: elapsed_max() may return too high a value on systems
  // where std::clock_t overflows or resets at surprising values.
  {
    return (double((std::numeric_limits<std::clock_t>::max)())
       - double(_start_time)) / double(CLOCKS_PER_SEC);
  }

  double elapsed_min() const            // return minimum value for elapsed()
   { return double(1)/double(CLOCKS_PER_SEC); }

 private:
  std::clock_t _start_time;
}; // timer

} // namespace test

#endif  // BOOST_PARAMETER_TEST_TIMER_HPP
