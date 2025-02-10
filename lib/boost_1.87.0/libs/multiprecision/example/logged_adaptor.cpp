///////////////////////////////////////////////////////////////
//  Copyright 2013 John Maddock. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at https://www.boost.org/LICENSE_1_0.txt

//[logged_adaptor

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

template <unsigned D>
inline void log_postfix_event(const mpfi_float_backend<D>& val, const char* event_description)
{
   // Print out the (relative) diameter of the interval:
   using namespace boost::multiprecision;
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
} // namespace boost::multiprecision
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

   for(unsigned i = 0; i < 13; ++i)
   {
      logged_type b = a * 9;
      b /= 10;
      a -= b;
   }
   std::cout << "Final value was: " << a << std::endl;
   return 0;
}

//]

/*
//[logged_adaptor_output

Diameter was -0 after operation: Default construct
Diameter was 0 after operation: Assignment from arithmetic type
Diameter was 3.34096e-51 after operation: /=
Diameter was -0 after operation: Default construct
Diameter was 5.93948e-51 after operation: *
Diameter was 7.42435e-51 after operation: /=
Diameter was 1.00229e-49 after operation: -=
Diameter was -0 after operation: Default construct
Diameter was 1.00229e-49 after operation: *
Diameter was 1.02085e-49 after operation: /=
Diameter was 1.92105e-48 after operation: -=
Diameter was -0 after operation: Default construct
Diameter was 1.92105e-48 after operation: *
Diameter was 1.92279e-48 after operation: /=
Diameter was 3.65156e-47 after operation: -=
Diameter was -0 after operation: Default construct
Diameter was 3.65156e-47 after operation: *
Diameter was 3.65163e-47 after operation: /=
Diameter was 6.93803e-46 after operation: -=
Diameter was -0 after operation: Default construct
Diameter was 6.93803e-46 after operation: *
Diameter was 6.93803e-46 after operation: /=
Diameter was 1.31823e-44 after operation: -=
Diameter was -0 after operation: Default construct
Diameter was 1.31823e-44 after operation: *
Diameter was 1.31823e-44 after operation: /=
Diameter was 2.50463e-43 after operation: -=
Diameter was -0 after operation: Default construct
Diameter was 2.50463e-43 after operation: *
Diameter was 2.50463e-43 after operation: /=
Diameter was 4.7588e-42 after operation: -=
Diameter was -0 after operation: Default construct
Diameter was 4.7588e-42 after operation: *
Diameter was 4.7588e-42 after operation: /=
Diameter was 9.04171e-41 after operation: -=
Diameter was -0 after operation: Default construct
Diameter was 9.04171e-41 after operation: *
Diameter was 9.04171e-41 after operation: /=
Diameter was 1.71793e-39 after operation: -=
Diameter was -0 after operation: Default construct
Diameter was 1.71793e-39 after operation: *
Diameter was 1.71793e-39 after operation: /=
Diameter was 3.26406e-38 after operation: -=
Diameter was -0 after operation: Default construct
Diameter was 3.26406e-38 after operation: *
Diameter was 3.26406e-38 after operation: /=
Diameter was 6.20171e-37 after operation: -=
Diameter was -0 after operation: Default construct
Diameter was 6.20171e-37 after operation: *
Diameter was 6.20171e-37 after operation: /=
Diameter was 1.17832e-35 after operation: -=
Diameter was -0 after operation: Default construct
Diameter was 1.17832e-35 after operation: *
Diameter was 1.17832e-35 after operation: /=
Diameter was 2.23882e-34 after operation: -=
Final value was: {1e-14,1e-14}
//]
*/
