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
namespace boost {
   namespace multiprecision {

      bool event_has_been_called = false;

      template <unsigned D>
      inline void log_postfix_event(const mpfi_float_backend<D>& val, const char* event_description)
      {
         // Print out the (relative) diameter of the interval:
         using namespace boost::multiprecision;
         event_has_been_called = true;
         number<mpfr_float_backend<D> > diam;
         mpfi_diam(diam.backend().data(), val.data());
         std::cout << "Diameter was " << diam << " after operation: " << event_description << std::endl;
      }
      template <unsigned D, class T>
      inline void log_postfix_event(const mpfi_float_backend<D>&, const T&, const char* event_description)
      {
         // This version is never called in this example.
      }

   }
}
//
// Now we can include the actual multiprecision headers and make the types concrete:
//
#include <boost/multiprecision/mpfi.hpp>
#include <boost/multiprecision/logged_adaptor.hpp>


int main()
{
   using namespace boost::multiprecision;
   using logged_type = logged_adaptor_t<mpfi_float_50>;
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

