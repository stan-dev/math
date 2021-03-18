#ifndef STAN_MATH_TEST_FIXTURE_ODE_SHO_HPP
#define STAN_MATH_TEST_FIXTURE_ODE_SHO_HPP

#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/rev/functor/coupled_mm.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

template <typename T>
struct harmonic_oscillator_ode_base {
  harm_osc_ode_fun f;
  harm_osc_ode_fun_eigen f_eigen;
  harm_osc_ode_data_fun f_data;
  harm_osc_ode_data_fun_eigen f_data_eigen;
  harm_osc_ode_wrong_size_1_fun f_wrong_size_1;
  harm_osc_ode_wrong_size_2_fun f_wrong_size_2;

  using T_t = std::tuple_element_t<2, T>;
  using T_y = std::tuple_element_t<3, T>;
  using T_p = std::tuple_element_t<4, T>;

  std::vector<T_p> theta;
  Eigen::Matrix<T_y, -1, 1> y0;
  double t0;
  std::vector<T_t> ts;
  std::vector<double> x_r;
  std::vector<int> x_i;
  double rtol;
  double atol;
  int max_num_step;

  harmonic_oscillator_ode_base()
      : theta{0.15},
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

  int dim() { return 2; }
  int param_size() { return 1; }
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
struct harmonic_oscillator_test
    : public harmonic_oscillator_ode_base<T>,
      public ODETestFixture<harmonic_oscillator_test<T>> {
  harmonic_oscillator_test() : harmonic_oscillator_ode_base<T>() {}

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

  void test_good() { ASSERT_NO_THROW(apply_solver()); }

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

    const double rtol_ = this->rtol;
    this->rtol = -1;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "relative_tolerance");
    this->rtol = rtol_;

    const double atol_ = this->atol;
    this->atol = -1;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "absolute_tolerance");
    this->atol = atol_;

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

/**
 * Bad ODE RHS type tests.
 *
 */
template <typename T>
struct harmonic_oscillator_bad_ode_test
    : public harmonic_oscillator_ode_base<T>,
      public ODETestFixture<harmonic_oscillator_bad_ode_test<T>> {
  harmonic_oscillator_bad_ode_test() : harmonic_oscillator_ode_base<T>() {}

  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(this->f_wrong_size_1, stan::math::to_array_1d(this->y0),
               this->t0, this->ts, this->theta, this->x_r, this->x_i, 0);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    return sol(this->f_wrong_size_1, stan::math::to_array_1d(this->y0),
               this->t0, this->ts, this->theta, this->x_r, this->x_i, 0,
               this->rtol, this->atol, this->max_num_step);
  }

  void test_bad_ode() {
    std::string error_msg = "dy_dt (3) and states (2) must match in size";
    EXPECT_THROW_MSG(apply_solver_tol(), std::invalid_argument, error_msg);
  }
};

/**
 * ODE RHS type that utilize <code>double</code> & <code>int</code> data
 *
 */
template <typename T>
struct harmonic_oscillator_data_test
    : public harmonic_oscillator_ode_base<T>,
      public ODETestFixture<harmonic_oscillator_data_test<T>> {
  harmonic_oscillator_data_test() : harmonic_oscillator_ode_base<T>() {
    this->x_r = std::vector<double>(3, 1);
    this->x_i = std::vector<int>(2, 0);
  }

  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(this->f_data_eigen, this->y0, this->t0, this->ts, 0, this->theta,
               this->x_r, this->x_i);
  }

  template <typename T1, typename T2>
  auto apply_solver(T1&& init, T2&& theta_in) {
    std::tuple_element_t<0, T> sol;
    return sol(this->f_data_eigen, init, this->t0, this->ts, nullptr, theta_in,
               this->x_r, this->x_i);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    return sol(this->f_data_eigen, this->y0, this->t0, this->ts, this->rtol,
               this->atol, this->max_num_step, 0, this->theta, this->x_r,
               this->x_i);
  }

  void test_bad_param_and_data() {
    const auto theta_ = this->theta;
    this->theta.resize(0);
    EXPECT_THROW_MSG(apply_solver(), std::out_of_range, "vector");
    this->theta = theta_;

    const auto x_r_ = this->x_r;
    this->x_r.resize(0);
    EXPECT_THROW_MSG(apply_solver(), std::out_of_range, "vector");
    this->x_r = x_r_;

    const auto x_i_ = this->x_i;
    this->x_i.resize(0);
    EXPECT_THROW_MSG(apply_solver(), std::out_of_range, "vector");
    this->x_i = x_i_;

    this->max_num_step = -1;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "max_num_steps");
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
