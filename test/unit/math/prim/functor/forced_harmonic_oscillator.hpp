#ifndef TEST_UNIT_MATH_ODE_FORCED_HARMONIC_OSCILLATOR
#define TEST_UNIT_MATH_ODE_FORCED_HARMONIC_OSCILLATOR

#include <stan/math/prim.hpp>
#include <stdexcept>
#include <vector>

struct forced_harm_osc_ode_fun {
  template <typename T0, typename T1, typename T2>
  inline std::vector<stan::return_type_t<T1, T2>>
  // initial time
  // initial positions
  // parameters
  // double data
  // integer data
  operator()(const T0& t_in, const std::vector<T1>& y_in,
             const std::vector<T2>& theta, const std::vector<double>& x,
             const std::vector<int>& x_int, std::ostream* msgs) const {
    if (y_in.size() != 2)
      throw std::domain_error(
          "this function was called with inconsistent state");

    std::vector<stan::return_type_t<T1, T2>> res;
    res.push_back(y_in.at(1));
    res.push_back(-y_in.at(0) * sin(theta.at(1) * t_in)
                  - theta.at(0) * y_in.at(1));

    return res;
  }
};

#endif
