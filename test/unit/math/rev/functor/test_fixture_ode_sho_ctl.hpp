#ifndef STAN_MATH_TEST_FIXTURE_ODE_SHO_CTL_HPP
#define STAN_MATH_TEST_FIXTURE_ODE_SHO_CTL_HPP

#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/rev/functor/coupled_mm.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode_sho.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

/**
 * Inheriting base type, various fixtures differs by the type of ODE functor
 * used in <code>apply_solver</code> calls, intended for different kind of
 * tests. This harmonic oscillator test class is intended for use with the
 * adjoint ODE solver which has additional control parameters.
 *
 */
template <typename T>
struct harmonic_oscillator_ctl_test
    : public harmonic_oscillator_ode_base<T>,
      public ODETestFixture<harmonic_oscillator_test<T>> {
  double rtol_f;
  double rtol_b;
  double rtol_q;
  Eigen::VectorXd atol_f;
  Eigen::VectorXd atol_b;
  double atol_q;
  int num_steps_check;
  int inter_poly;
  int solv_f;
  int solv_b;

  harmonic_oscillator_ctl_test()
      : harmonic_oscillator_ode_base<T>(),
        rtol_f(this->rtol / 8.0),
        rtol_b(this->rtol / 4.0),
        rtol_q(this->rtol),
        atol_f(Eigen::VectorXd::Constant(this->y0.size(), this->atol / 6.0)),
        atol_b(Eigen::VectorXd::Constant(this->y0.size(), this->atol / 3.0)),
        atol_q(this->atol),
        num_steps_check(100),
        inter_poly(CV_HERMITE),
        solv_f(CV_BDF),
        solv_b(CV_ADAMS) {}

  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(this->f_eigen, this->y0, this->t0, this->ts, nullptr,
               this->theta, this->x_r, this->x_i);
  }

  template <typename T1, typename T2>
  auto apply_solver(T1&& init, T2&& theta_in) {
    std::tuple_element_t<0, T> sol;
    return sol(this->f_eigen, init, this->t0, this->ts, nullptr, theta_in,
               this->x_r, this->x_i);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    return sol(this->f_eigen, this->y0, this->t0, this->ts, this->rtol,
               this->atol, this->max_num_step, nullptr, this->theta, this->x_r,
               this->x_i);
  }

  auto apply_solver_tol_ctl() {
    std::tuple_element_t<0, T> sol;
    return sol(this->f_eigen, this->y0, this->t0, this->ts, this->rtol_f,
               this->atol_f, this->rtol_b, this->atol_b, this->rtol_q,
               this->atol_q, this->max_num_step, this->num_steps_check,
               this->inter_poly, this->solv_f, this->solv_b, nullptr,
               this->theta, this->x_r, this->x_i);
  }

  void test_bad() {
    const auto y0_(this->y0);

    this->y0.resize(0);
    EXPECT_THROW_MSG(apply_solver(), std::invalid_argument,
                     "initial state has size 0");
    this->y0 = y0_;

    const auto t0_ = this->t0;
    this->t0 = 2.0;
    EXPECT_THROW_MSG(apply_solver(), std::domain_error,
                     "initial time is 2, but must be less than 0.1");
    this->t0 = t0_;

    const auto ts_ = this->ts;
    this->ts.resize(0);
    EXPECT_THROW_MSG(apply_solver(), std::invalid_argument, "times has size 0");
    this->ts = ts_;

    this->ts.resize(2);
    this->ts[0] = 3.0;
    this->ts[1] = 1.0;
    EXPECT_THROW_MSG(apply_solver(), std::domain_error,
                     "times is not a valid sorted vector");
    this->ts = ts_;

    const double rtol_f_ = this->rtol_f;
    this->rtol_f = -1;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::domain_error,
                     "relative_tolerance_forward");
    this->rtol_f = rtol_f_;

    const double atol_f_ = this->atol_f(0);
    this->atol_f(0) = -1;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::domain_error,
                     "absolute_tolerance_forward");
    this->atol_f(0) = atol_f_;

    this->atol_f.resize(1);
    this->atol_f(0) = atol_f_;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::invalid_argument,
                     "absolute_tolerance_forward");
    this->atol_f.resize(2);
    this->atol_f(0) = atol_f_;
    this->atol_f(1) = atol_f_;

    const double rtol_b_ = this->rtol_b;
    this->rtol_b = -1;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::domain_error,
                     "relative_tolerance_backward");
    this->rtol_b = rtol_b_;

    const double atol_b_ = this->atol_b(0);
    this->atol_b(0) = -1;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::domain_error,
                     "absolute_tolerance_backward");
    this->atol_b(0) = atol_b_;

    this->atol_b.resize(1);
    this->atol_b(0) = atol_b_;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::invalid_argument,
                     "absolute_tolerance_backward");
    this->atol_b.resize(2);
    this->atol_b(0) = atol_b_;
    this->atol_b(1) = atol_b_;

    const double rtol_q_ = this->rtol_q;
    this->rtol_q = -1;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::domain_error,
                     "relative_tolerance_quadrature");
    this->rtol_q = rtol_q_;

    const double atol_q_ = this->atol_q;
    this->atol_q = -1;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::domain_error,
                     "absolute_tolerance_quadrature");
    this->atol_q = atol_q_;

    const int max_num_step_ = this->max_num_step;
    this->max_num_step = -1;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::domain_error,
                     "max_num_steps");
    this->max_num_step = max_num_step_;

    const int num_steps_check_ = this->num_steps_check;
    this->num_steps_check = -1;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::domain_error,
                     "num_steps_between_checkpoints");
    this->num_steps_check = num_steps_check_;

    const int inter_poly_ = this->inter_poly;
    this->inter_poly = 0;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::invalid_argument,
                     "interpolation_polynomial");
    this->inter_poly = inter_poly_;

    const int solv_f_ = this->solv_f;
    this->solv_f = 0;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::invalid_argument,
                     "solver_forward");
    this->solv_f = solv_f_;

    const int solv_b_ = this->solv_b;
    this->solv_b = 0;
    EXPECT_THROW_MSG(apply_solver_tol_ctl(), std::invalid_argument,
                     "solver_backward");
    this->solv_b = solv_b_;

    const auto theta_ = this->theta;
    const auto x_r_ = this->x_r;
    const auto x_i_ = this->x_i;

    // NaN errors
    double nan = std::numeric_limits<double>::quiet_NaN();
    std::stringstream expected_is_nan;
    expected_is_nan << "is " << nan;

    this->y0[0] = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial state");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_nan.str());
    this->y0 = y0_;

    this->t0 = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial time");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_nan.str());
    this->t0 = t0_;

    this->ts[0] = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "times");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_nan.str());
    this->ts = ts_;

    this->theta[0] = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "ode parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_nan.str());
    this->theta = theta_;

    this->x_r.push_back(nan);
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "ode parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_nan.str());
    this->x_r = x_r_;

    // inf test
    std::stringstream expected_is_inf;
    expected_is_inf << "is " << std::numeric_limits<double>::infinity();
    std::stringstream expected_is_neg_inf;
    expected_is_neg_inf << "is " << -std::numeric_limits<double>::infinity();
    double inf = std::numeric_limits<double>::infinity();

    this->y0[0] = inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial state");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_inf.str());
    this->y0[0] = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial state");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_neg_inf.str());
    this->y0 = y0_;

    this->t0 = inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial time");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_inf.str());
    this->t0 = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial time");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_neg_inf.str());
    this->t0 = t0_;

    this->ts.back() = inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "times");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_inf.str());
    this->ts.back() = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "times");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_neg_inf.str());
    this->ts = ts_;

    this->theta[0] = inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "ode parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_inf.str());
    this->theta[0] = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "ode parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_neg_inf.str());
    this->theta = theta_;

    this->x_r = std::vector<double>{inf};
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "ode parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_inf.str());
    this->x_r[0] = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "ode parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_neg_inf.str());
    this->x_r = x_r_;
  }

  void test_value(double t0_in) {
    this->t0 = t0_in;
    for (size_t i = 0; i < this->ts.size(); ++i) {
      this->ts[i] = this->t0 + 0.1 * (i + 1);
    }

    this->rtol = 1e-8;
    this->atol = 1e-10;
    this->max_num_step = 1e6;
    auto res = apply_solver_tol();

    EXPECT_NEAR(0.995029, stan::math::value_of(res[0][0]), 1e-5);
    EXPECT_NEAR(-0.0990884, stan::math::value_of(res[0][1]), 1e-5);

    EXPECT_NEAR(-0.421907, stan::math::value_of(res[99][0]), 1e-5);
    EXPECT_NEAR(0.246407, stan::math::value_of(res[99][1]), 1e-5);
  }
};

#endif
