#ifndef STAN_MATH_TEST_FIXTURE_DAE_LINEAR_HPP
#define STAN_MATH_TEST_FIXTURE_DAE_LINEAR_HPP

#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/coupled_mm.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

/** 
 * Linear constant coef DAE in the form
 * x''(t) = z(t)
 * x(t) + z(t) = 3 sin(t)
 *
 * with IC x(pi) = 1.0, x'(pi) = 0.0
 *
 */
template <typename T>
struct linear_dae_base {
  struct func {
  template <typename T0, typename Tyy, typename Typ, typename Tpar>
  inline Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> operator()(
      const T0& t, const Eigen::Matrix<Tyy, -1, 1>& yy, const
      Eigen::Matrix<Typ, -1, 1> & yp,
      std::ostream* msgs, const Tpar& theta) const {
    if (yy.size() != 3 || yp.size() != 3)
      throw std::domain_error(
          "this function was called with inconsistent state");

    Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> res(3);

    auto yy1 = yy(0);
    auto yy2 = yy(1);
    auto yy3 = yy(2);

    auto yp1 = yp(0);
    auto yp2 = yp(1);

    res[0] = yp1 - yy2;
    res[1] = yp2 - yy3;
    res[2] = yy1 + yy3 - theta * sin(t);

    return res;
  }
};

  func f;

  using T_t = std::tuple_element_t<2, T>;
  using Tyy = std::tuple_element_t<3, T>;
  using Typ = std::tuple_element_t<4, T>;
  using T_p = std::tuple_element_t<5, T>;

  T_p theta;
  Eigen::Matrix<Tyy, -1, 1> yy0;
  Eigen::Matrix<Typ, -1, 1> yp0;
  double t0;
  std::vector<T_t> ts;
  double rtol;
  double atol;
  int max_num_step;

  linear_dae_base()
    : theta(3.0),
      yy0(3),
      yp0(3),
      t0(3.14159265358979323846),
      ts{2.4 * t0},
      rtol(1.e-5),
      atol(1.e-8),
      max_num_step(1000) {
    yy0 << 1.0, 0.0, -1.0;
    yp0 << 0.0, -1.0, 0.0;
  }

  std::vector<T_t>& times() { return ts; }
  Eigen::VectorXd init() {
    Eigen::VectorXd joined_init(yy0.size() + yp0.size());
    joined_init << stan::math::value_of(yy0), stan::math::value_of(yp0);
    return joined_init;
  }
  std::vector<double> param() { return {theta}; }
};

template <typename T>
struct linear_dae_functor : public linear_dae_base<T> {
  template <typename Tx>
  Eigen::Matrix<Tx, -1, 1> operator()(const Eigen::Matrix<Tx, -1, 1>& x) const {
    std::tuple_element_t<0, T> sol;
    auto ys = sol(this->f, this->yy0, this->yp0, this->t0, this->ts,
                  nullptr, x[0]);
    return ys[0];
  }
};

template <typename T>
struct linear_dae_test
    : public linear_dae_base<T>,
      public ODETestFixture<linear_dae_test<T>> {
  linear_dae_functor<T> dae_sol;

  /** 
   * analytical solution based on theta = 3.0
   * 
   */
  auto analy_sol_functor() {
    auto f = [&](double t) {
      const double pi = 3.14159265358979323846;
      Eigen::VectorXd y(3);
      y << -cos(t) + 1.5 * pi * cos(t) - 1.5 * t * cos(t) + 1.5 * sin(t),
        sin(t) - 1.5 * pi * sin(t) - 1.5 * cos(t) + 1.5 * t * sin(t) +
        1.5 * cos(t),
        cos(t) - 1.5 * pi * cos(t) + 1.5 * t * cos(t) + 1.5 * sin(t);
      return y;
    };
    return f;
  }

  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(this->f, this->yy0, this->yp0, this->t0, this->ts, nullptr,
               this->theta);
  }

  template <typename T1, typename T2>
  auto apply_solver(T1&& init, T2&& theta_in) {
    const int n = init.size()/2;
    std::tuple_element_t<0, T> sol;
    return sol(this->f, init.head(n), init.tail(n), this->t0,
               this->ts, nullptr, theta_in[0]);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    return sol(this->f, this->yy0, this->yp0, this->t0, this->ts, this->rtol,
               this->atol, this->max_num_step, nullptr, this->theta);
  }
};

#endif
