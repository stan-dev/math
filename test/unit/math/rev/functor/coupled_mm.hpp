#ifndef TEST_UNIT_MATH_REV_ARR_FUNCTOR_COUPLED_MM_HPP
#define TEST_UNIT_MATH_REV_ARR_FUNCTOR_COUPLED_MM_HPP

#include <stan/math/rev/core.hpp>
#include <vector>

struct coupled_mm_ode_fun {
  template <typename T0, typename T1, typename T2>
  inline std::vector<typename stan::return_type<T1, T2>::type>
  // initial time
  // initial positions
  // parameters
  // double data
  // integer data
  operator()(const T0& t_in, const std::vector<T1>& y,
             const std::vector<T2>& parms, const std::vector<double>& sx,
             const std::vector<int>& sx_int, std::ostream* msgs) const {
    std::vector<typename stan::return_type<T1, T2>::type> ydot(2);

    const T2 act = parms[0];
    const T2 KmA = parms[1];
    const T2 deact = parms[2];
    const T2 KmAp = parms[3];

    ydot[0]
        = -1 * (act * y[0] / (KmA + y[0])) + 1 * (deact * y[1] / (KmAp + y[1]));
    ydot[1]
        = 1 * (act * y[0] / (KmA + y[0])) - 1 * (deact * y[1] / (KmAp + y[1]));

    return (ydot);
  }
};

#endif
