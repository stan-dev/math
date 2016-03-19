#pragma once

#include <stan/math/prim/arr/functor/integrate_ode.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/scal/fun/pow.hpp>

// defines taken over from cOde generated code
#define Vm1 parms[0]
 #define Km1 parms[1]
 #define Vm2 parms[2]
 #define Km2 parms[3]
 #define k3 parms[4]
 #define Km3 parms[5]
 #define Vm4 parms[6]
 #define Km4 parms[7]
 #define k5 parms[8]
 #define Km5 parms[9]
 #define Vm6 parms[10]
 #define Km6 parms[11]
 #define k7 parms[12]
 #define Km7 parms[13]
 #define Vm8 parms[14]
 #define Km8 parms[15]
 #define Inh parms[16]
 #define Ki8 parms[17]


struct hornberg_ode_fun {
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
    std::vector<typename stan::return_type<T1,T2>::type> ydot(8);

    ydot[0] = -1*((Vm1*y[0])/(Km1+y[0]))+1*((Vm2*y[1])/(Km2+y[1]));
    ydot[1] = 1*((Vm1*y[0])/(Km1+y[0]))-1*((Vm2*y[1])/(Km2+y[1]));
    ydot[2] = -1*((k3*y[0]*y[2])/(Km3+y[2]))+1*((Vm4*y[3])/(Km4+y[3]));
    ydot[3] = 1*((k3*y[0]*y[2])/(Km3+y[2]))-1*((Vm4*y[3])/(Km4+y[3]));
    ydot[4] = -1*((k5*y[3]*y[4])/(Km5+y[4]))+1*((Vm6*y[5])/(Km6+y[5]));
    ydot[5] = 1*((k5*y[3]*y[4])/(Km5+y[4]))-1*((Vm6*y[5])/(Km6+y[5]));
    ydot[6] = -1*((k7*y[5]*y[6])/(Km7+y[6]))+1*((Vm8*y[7])/(Km8*(1+Inh/Ki8+y[7]/Km8)));
    ydot[7] = 1*((k7*y[5]*y[6])/(Km7+y[6]))-1*((Vm8*y[7])/(Km8*(1+Inh/Ki8+y[7]/Km8)));

    return(ydot);
  }
};
