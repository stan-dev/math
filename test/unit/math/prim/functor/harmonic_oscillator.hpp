#ifndef TEST_UNIT_MATH_ODE_HARMONIC_OSCILLATOR
#define TEST_UNIT_MATH_ODE_HARMONIC_OSCILLATOR

#include <stan/math/prim.hpp>
#include <stdexcept>
#include <vector>

struct harm_osc_ode_fun {
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
    res.push_back(-y_in.at(0) - theta.at(0) * y_in.at(1));

    return res;
  }
};

struct harm_osc_ode_fun_eigen {
  template <typename T0, typename T1, typename T2>
  inline auto operator()(const T0& t_in, const T1& y_in, std::ostream* msgs,
                         const std::vector<T2>& theta,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int) const {
    if (y_in.size() != 2)
      throw std::domain_error(
          "this function was called with inconsistent state");

    Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1> res(2);
    res(0) = y_in(1);
    res(1) = -y_in(0) - theta[0] * y_in(1);

    return res;
  }
};

struct harm_osc_ode_data_fun {
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
    res.push_back(x.at(0) * y_in.at(1) + x_int.at(0));
    res.push_back(-x.at(1) * y_in.at(0) - x.at(2) * theta.at(0) * y_in.at(1)
                  + x_int.at(1));

    return res;
  }
};

struct harm_osc_ode_data_fun_eigen {
  template <typename T0, typename T1, typename T2>
  inline auto operator()(const T0& t_in, const T1& y_in, std::ostream* msgs,
                         const std::vector<T2>& theta,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int) const {
    if (y_in.size() != 2)
      throw std::domain_error(
          "this function was called with inconsistent state");

    const T2& p = theta.at(0);

    Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> res(2);
    res << (x.at(0) * y_in(1) + x_int.at(0)),
        (-x.at(1) * y_in(0) - x.at(2) * p * y_in(1) + x_int.at(1));

    return res;
  }
};

struct harm_osc_ode_wrong_size_1_fun {
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
    res.push_back(-y_in.at(0) - theta.at(0) * y_in.at(1));
    res.push_back(0);

    return res;
  }
};

struct harm_osc_ode_wrong_size_2_fun {
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
    res.push_back(0);

    return res;
  }
};

#endif
