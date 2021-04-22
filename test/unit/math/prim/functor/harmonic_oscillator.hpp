#ifndef TEST_UNIT_MATH_ODE_HARMONIC_OSCILLATOR
#define TEST_UNIT_MATH_ODE_HARMONIC_OSCILLATOR

#include <stan/math/prim.hpp>
#include <stdexcept>
#include <vector>

struct harm_osc_ode_fun {
  // initial time
  // initial positions
  // parameters
  // double data
  // integer data
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  inline std::vector<stan::return_type_t<T1, T2>>
  operator()(const T0& t_in, const T1& y_in,
             const T2& theta, const T3& x,
             const T4& x_int, std::ostream* msgs) const {
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
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  inline auto operator()(const T0& t_in,
                         const T1& y_in,
                         std::ostream* msgs, const T2& theta,
                         const T3& x,
                         const T4& x_int) const {
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
  // initial time
  // initial positions
  // parameters
  // double data
  // integer data
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  inline std::vector<stan::return_type_t<T1, T2>>
  operator()(const T0& t_in, const T1& y_in,
             const T2& theta, const T3& x,
             const T4& x_int, std::ostream* msgs) const {
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
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  inline auto operator()(const T0& t_in,
                         const T1& y_in,
                         std::ostream* msgs, const T2& theta,
                         const T3& x,
                         const T4& x_int) const {
    if (y_in.size() != 2)
      throw std::domain_error(
          "this function was called with inconsistent state");

    const auto& p = theta.at(0);

    Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> res(2);
    res << (x.at(0) * y_in(1) + x_int.at(0)),
        (-x.at(1) * y_in(0) - x.at(2) * p * y_in(1) + x_int.at(1));

    return res;
  }
};

struct harm_osc_ode_wrong_size_1_fun {
  // initial time
  // initial positions
  // parameters
  // double data
  // integer data
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  inline std::vector<stan::return_type_t<T1, T2>>
  operator()(const T0& t_in, const T1& y_in,
             const T2& theta, const T3& x,
             const T4& x_int, std::ostream* msgs) const {
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
  // initial time
  // initial positions
  // parameters
  // double data
  // integer data
  template <typename T0, typename T1, typename T2, typename T3, typename T4>
  inline std::vector<stan::return_type_t<T1, T2>>
  operator()(const T0& t_in, const T1& y_in,
             const T2& theta, const T3& x,
             const T4& x_int, std::ostream* msgs) const {
    if (y_in.size() != 2)
      throw std::domain_error(
          "this function was called with inconsistent state");

    std::vector<stan::return_type_t<T1, T2>> res;
    res.push_back(0);

    return res;
  }
};

#endif
