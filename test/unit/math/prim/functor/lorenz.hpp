#ifndef TEST__UNIT__MATH__ODE__LORENZ_HPP
#define TEST__UNIT__MATH__ODE__LORENZ_HPP

#include <vector>

template <typename T0, typename T1, typename T2>
inline std::vector<stan::return_type_t<T1, T2>>
// initial time
// initial positions
// parameters
// double data
// integer data
lorenz_ode(const T0& t_in, const std::vector<T1>& y_in,
           const std::vector<T2>& theta, const std::vector<double>& x,
           const std::vector<int>& x_int) {
  std::vector<stan::return_type_t<T1, T2>> res;
  res.push_back(theta.at(0) * (y_in.at(1) - y_in.at(0)));
  res.push_back(theta.at(1) * y_in.at(0) - y_in.at(1)
                - y_in.at(0) * y_in.at(2));
  res.push_back(-theta.at(2) * y_in.at(2) + y_in.at(0) * y_in.at(1));
  return res;
}

struct lorenz_ode_fun {
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
    return lorenz_ode(t_in, y_in, theta, x, x_int);
  }
};

#endif
