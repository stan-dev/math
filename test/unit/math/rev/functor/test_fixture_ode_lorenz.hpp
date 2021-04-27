#ifndef STAN_MATH_TEST_FIXTURE_ODE_LORENZ_HPP
#define STAN_MATH_TEST_FIXTURE_ODE_LORENZ_HPP

#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/functor/lorenz.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

struct lorenz_ode_base {
  struct lorenz_rhs {
    template <typename T0, typename T1, typename T2>
    inline Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> operator()(
        const T0& t_in, const T1& y_in, std::ostream* msgs,
        const T2& theta) const {
      Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> res(3);
      res << theta.at(0) * (y_in(1) - y_in(0)),
          theta.at(1) * y_in(0) - y_in(1) - y_in(0) * y_in(2),
          -theta.at(2) * y_in(2) + y_in(0) * y_in(1);
      return res;
    }
  };

  lorenz_rhs f;
  Eigen::VectorXd y0;
  std::vector<double> theta;
  double t0;
  std::vector<double> ts;
  double rtol;
  double atol;
  int max_num_step;

  lorenz_ode_base()
      : f(),
        y0(3),
        theta{10.0, 28.0, 8.0 / 3.0},
        t0(0.0),
        ts(100),
        rtol(1.e-8),
        atol(1.e-10),
        max_num_step(100000) {
    y0 << 10.0, 1.0, 1.0;
    for (int i = 0; i < 100; i++)
      ts[i] = t0 + (0.1 * (i + 1));
  }

  std::vector<double>& times() { return ts; }
  Eigen::VectorXd init() { return y0; }
  std::vector<double> param() { return theta; }
};

/**
 * Inheriting base type, various fixtures differs by the type of ODE
 * functor used in <code>apply_solver</code> calls, intended for
 * different kind of tests.
 *
 */
template <typename T>
struct lorenz_test : public lorenz_ode_base,
                     public ODETestFixture<lorenz_test<T>> {
  lorenz_test() : lorenz_ode_base() {}

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
};

#endif
