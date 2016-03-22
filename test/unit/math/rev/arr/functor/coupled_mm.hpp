#pragma once

#include <stan/math/prim/arr/functor/integrate_ode.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/scal/fun/pow.hpp>

#define act parms[0] 
 #define KmA parms[1] 
 #define deact parms[2] 
 #define KmAp parms[3] 

struct coupled_mm_ode_fun {
  template <typename T0, typename T1, typename T2>
  inline
  std::vector<typename stan::return_type<T1,T2>::type>
  operator()(const T0& t_in, // initial time
             const std::vector<T1>& y, //initial positions
             const std::vector<T2>& parms, // parameters
             const std::vector<double>& sx, // double data
             const std::vector<int>& sx_int,
             std::ostream* msgs) const { // integer data
    using std::pow;
    std::vector<typename stan::return_type<T1,T2>::type> ydot(2);

    ydot[0] = -1*(act*y[0]/(KmA+y[0]))+1*(deact*y[1]/(KmAp+y[1]));
    ydot[1] = 1*(act*y[0]/(KmA+y[0]))-1*(deact*y[1]/(KmAp+y[1]));

    return(ydot);
  }
};
