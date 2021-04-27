#ifndef STAN_MATH_TEST_FIXTURE_ODE_FHO_HPP
#define STAN_MATH_TEST_FIXTURE_ODE_FHO_HPP

#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

template <typename T>
struct forced_harm_osc_base {
  struct forced_harm_osc {
    template <typename T0, typename T1, typename T2>
    inline Eigen::Matrix<stan::return_type_t<T0, T1, T2>, -1, 1> operator()(
        const T0& t_in, const T1& y_in, std::ostream* msgs,
        const T2& theta) const {
      if (y_in.size() != 2)
        throw std::domain_error(
            "this function was called with inconsistent state");

      Eigen::Matrix<stan::return_type_t<T0, T1, T2>, -1, 1> res(2);
      res << y_in(1),
          -y_in(0) * sin(theta.at(1) * t_in) - theta.at(0) * y_in(1);

      return res;
    }
  };

  forced_harm_osc f;

  using T_t = std::tuple_element_t<2, T>;
  using T_y = std::tuple_element_t<3, T>;
  using T_p = std::tuple_element_t<4, T>;

  std::vector<T_p> theta;
  Eigen::Matrix<T_y, -1, 1> y0;
  double t0;
  std::vector<T_t> ts;
  double rtol;
  double atol;
  int max_num_step;

  forced_harm_osc_base()
      : theta{0.15, 0.25},
        y0(2),
        t0(0),
        ts(100),
        rtol(1.e-6),
        atol(1.e-8),
        max_num_step(100000) {
    y0 << 1.0, 0.0;
    for (size_t i = 0; i < ts.size(); ++i) {
      ts[i] = t0 + 0.1 * (i + 1);
    }
  }

  std::vector<T_t>& times() { return ts; }
  Eigen::Matrix<T_y, -1, 1> init() { return y0; }
  std::vector<T_p> param() { return theta; }
};

/**
 * Inheriting base type, various fixtures differs by the type of ODE
 * functor used in <code>apply_solver</code> calls, intended for
 * different kind of tests.
 *
 */
template <typename T>
struct forced_harm_osc_ts_test
    : public forced_harm_osc_base<T>,
      public ODETestFixture<forced_harm_osc_ts_test<T>> {
  forced_harm_osc_ts_test() : forced_harm_osc_base<T>() {}

  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(this->f, this->y0, this->t0, this->ts, nullptr, this->theta);
  }

  template <typename T1, typename T2>
  auto apply_solver(T1&& init, T2&& theta_in) {
    std::tuple_element_t<0, T> sol;
    return sol(this->f, init, this->t0, this->ts, nullptr, theta_in);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    return sol(this->f, this->y0, this->t0, this->ts, this->rtol, this->atol,
               this->max_num_step, nullptr, this->theta);
  }

  template <typename Ttime, typename T1>
  auto eval_rhs(const Ttime& t, const T1& y) {
    return this->f(t, y, nullptr, this->theta);
  }
};

#endif
