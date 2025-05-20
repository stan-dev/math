///////////////////////////////////////////////////////////////
//  Copyright 2023 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt

#include <iostream>
#include <iomanip>
#include <boost/multiprecision/fwd.hpp>
//
// Begin by overloading log_postfix_event so we can capture each arithmetic event as it happens,
// unfortunately this must occur BEFORE we include the full header, so just include the forward
// declarations and define our overloads for now.  Note that in some cases we may need to just
// declare the overloads here, and define them once the types become concrete:
//
namespace boost::multiprecision {
      void log_postfix_event(const cpp_bin_float<113, backends::digit_base_2, void, std::int16_t, -16382, 16383>& val, const char* event_description);
}
//
// Now we can include the actual multiprecision headers and make the types concrete:
//
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/logged_adaptor.hpp>

namespace boost::multiprecision {

      bool event_has_been_called = false;

      inline void log_postfix_event(const cpp_bin_float<113, backends::digit_base_2, void, std::int16_t, -16382, 16383>& val, const char* event_description)
      {
         // Print out the (relative) diameter of the interval:
         using namespace boost::multiprecision;
         event_has_been_called = true;
      }
}

int main()
{
   using namespace boost::multiprecision;
   using logged_type = logged_adaptor_t<cpp_bin_float_quad>;
   //
   // Test case deliberately introduces cancellation error, relative size of interval
   // gradually gets larger after each operation:
   //
   logged_type a = 1;
   a /= 10;

   for (unsigned i = 0; i < 13; ++i)
   {
      logged_type b = a * 9;
      b /= 10;
      a -= b;
   }
   std::cout << "Final value was: " << a << std::endl;
   return boost::multiprecision::event_has_been_called ? 0 : -1;
}

