#ifndef TEST_UNIT_MATH_ODE_FN_OSCILLATOR
#define TEST_UNIT_MATH_ODE_FN_OSCILLATOR

#include <stan/math/prim.hpp>
#include <stdexcept>
#include <vector>

struct Field_Noyes_osc {
  template <typename T0, typename T1, typename T2>
  inline auto operator()(const T0& t_in,
                         const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
                         std::ostream* msgs, const std::vector<T2>& theta,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int) const {
    if (y.size() != 3)
      throw std::domain_error(
          "this function was called with inconsistent state");

    const T2& a = theta[0];

    Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> res(3);
    res(0) = a * (y(1) - y(0) * y(1) + y(0) - 8.375e6 * y(0) * y(0));
    res(1) = (-y(1) - y(0) * y(1) + y(2)) / a;
    res(2) = 0.161 * (y(0) - y(2));

    return res;
  }
};

#endif
