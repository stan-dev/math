#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/functor/lorenz.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>

std::ostream* msgs = nullptr;

template <typename F, typename... Args>
void sho_death_test(F harm_osc, std::vector<double>& y0, double t0,
                    std::vector<double>& ts, Args&... args) {
  EXPECT_DEATH(stan::math::integrate_ode_rk45(harm_osc, y0, t0, ts, args...),
               "");
}

template <typename F, typename... Args>
void sho_value_test(F harm_osc, std::vector<double>& y0, double t0,
                    std::vector<double>& ts, Args&... args) {
  std::vector<std::vector<double>> ode_res_vd
      = stan::math::integrate_ode_rk45(harm_osc, y0, t0, ts, args...);
  EXPECT_NEAR(0.995029, ode_res_vd[0][0], 1e-5);
  EXPECT_NEAR(-0.0990884, ode_res_vd[0][1], 1e-5);

  EXPECT_NEAR(-0.421907, ode_res_vd[99][0], 1e-5);
  EXPECT_NEAR(0.246407, ode_res_vd[99][1], 1e-5);
}

void sho_test(double t0) {
  harm_osc_ode_fun harm_osc;

  std::vector<double> theta;
  theta.push_back(0.15);

  std::vector<double> y0;
  y0.push_back(1.0);
  y0.push_back(0.0);

  std::vector<double> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x;
  std::vector<int> x_int;

  sho_value_test(harm_osc, y0, t0, ts, theta, x, x_int);
}

void sho_data_test(double t0) {
  harm_osc_ode_data_fun harm_osc;

  std::vector<double> theta;
  theta.push_back(0.15);

  std::vector<double> y0;
  y0.push_back(1.0);
  y0.push_back(0.0);

  std::vector<double> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  sho_value_test(harm_osc, y0, t0, ts, theta, x, x_int);
}

TEST(StanMathOde_integrate_ode_rk45, harmonic_oscillator) {
  sho_test(0.0);
  sho_test(1.0);
  sho_test(-1.0);

  sho_data_test(0.0);
  sho_data_test(1.0);
  sho_data_test(-1.0);
}

TEST(StanMathOde_integrate_ode_rk45, error_conditions) {
  using stan::math::integrate_ode_rk45;
  harm_osc_ode_data_fun harm_osc;

  std::vector<double> theta;
  theta.push_back(0.15);

  std::vector<double> y0;
  y0.push_back(1.0);
  y0.push_back(0.0);

  double t0 = 0;

  std::vector<double> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  ASSERT_NO_THROW((integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x, x_int)));

  std::vector<double> y0_bad;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0_bad, t0, ts, theta, x, x_int)),
      std::invalid_argument, "initial state has size 0");

  double t0_bad = 2.0;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0_bad, ts, theta, x, x_int)),
      std::domain_error, "initial time is 2, but must be less than 0.1");

  std::vector<double> ts_bad;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts_bad, theta, x, x_int)),
      std::invalid_argument, "times has size 0");

  ts_bad.push_back(3);
  ts_bad.push_back(1);
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts_bad, theta, x, x_int)),
      std::domain_error, "times is not a valid sorted vector");

  std::vector<double> theta_bad;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts, theta_bad, x, x_int)),
      std::out_of_range, "vector");

  std::vector<double> x_bad;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x_bad, x_int)),
      std::out_of_range, "vector");

  std::vector<int> x_int_bad;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x, x_int_bad)),
      std::out_of_range, "vector");

  EXPECT_THROW_MSG((integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x, x_int, 0,
                                       -1, 1e-6, 10)),
                   std::domain_error, "relative_tolerance");

  EXPECT_THROW_MSG((integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x, x_int, 0,
                                       1e-6, -1, 10)),
                   std::domain_error, "absolute_tolerance");

  EXPECT_THROW_MSG((integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x, x_int, 0,
                                       1e-6, 1e-6, -1)),
                   std::domain_error, "max_num_steps");
}

TEST(StanMathOde_integrate_ode_rk45, error_conditions_nan) {
  using stan::math::integrate_ode_rk45;
  harm_osc_ode_data_fun harm_osc;

  std::vector<double> theta;
  theta.push_back(0.15);

  std::vector<double> y0;
  y0.push_back(1.0);
  y0.push_back(0.0);

  double t0 = 0;

  std::vector<double> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  ASSERT_NO_THROW((integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x, x_int)));

  double nan = std::numeric_limits<double>::quiet_NaN();
  std::stringstream expected_is_nan;
  expected_is_nan << "is " << nan;

  std::vector<double> y0_bad = y0;
  y0_bad[0] = nan;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0_bad, t0, ts, theta, x, x_int)),
      std::domain_error, "initial state");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0_bad, t0, ts, theta, x, x_int)),
      std::domain_error, expected_is_nan.str());

  double t0_bad = nan;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0_bad, ts, theta, x, x_int)),
      std::domain_error, "initial time");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0_bad, ts, theta, x, x_int)),
      std::domain_error, expected_is_nan.str());

  std::vector<double> ts_bad = ts;
  ts_bad[0] = nan;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts_bad, theta, x, x_int)),
      std::domain_error, "times");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts_bad, theta, x, x_int)),
      std::domain_error, expected_is_nan.str());

  std::vector<double> theta_bad = theta;
  theta_bad[0] = nan;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts, theta_bad, x, x_int)),
      std::domain_error, "parameters and data");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts, theta_bad, x, x_int)),
      std::domain_error, expected_is_nan.str());

  if (x.size() > 0) {
    std::vector<double> x_bad = x;
    x_bad[0] = nan;
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x_bad, x_int)),
        std::domain_error, "parameters and data");
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x_bad, x_int)),
        std::domain_error, expected_is_nan.str());
  }
}

TEST(StanMathOde_integrate_ode_rk45, error_conditions_inf) {
  using stan::math::integrate_ode_rk45;
  harm_osc_ode_data_fun harm_osc;
  std::stringstream expected_is_inf;
  expected_is_inf << "is " << std::numeric_limits<double>::infinity();
  std::stringstream expected_is_neg_inf;
  expected_is_neg_inf << "is " << -std::numeric_limits<double>::infinity();

  std::vector<double> theta;
  theta.push_back(0.15);

  std::vector<double> y0;
  y0.push_back(1.0);
  y0.push_back(0.0);

  double t0 = 0;

  std::vector<double> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  ASSERT_NO_THROW((integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x, x_int)));

  double inf = std::numeric_limits<double>::infinity();
  std::vector<double> y0_bad = y0;
  y0_bad[0] = inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0_bad, t0, ts, theta, x, x_int)),
      std::domain_error, "initial state");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0_bad, t0, ts, theta, x, x_int)),
      std::domain_error, expected_is_inf.str());

  y0_bad[0] = -inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0_bad, t0, ts, theta, x, x_int)),
      std::domain_error, "initial state");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0_bad, t0, ts, theta, x, x_int)),
      std::domain_error, expected_is_neg_inf.str());

  double t0_bad = inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0_bad, ts, theta, x, x_int)),
      std::domain_error, "initial time");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0_bad, ts, theta, x, x_int)),
      std::domain_error, expected_is_inf.str());
  t0_bad = -inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0_bad, ts, theta, x, x_int)),
      std::domain_error, "initial time");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0_bad, ts, theta, x, x_int)),
      std::domain_error, expected_is_neg_inf.str());

  std::vector<double> ts_bad = ts;
  ts_bad.push_back(inf);
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts_bad, theta, x, x_int)),
      std::domain_error, "times");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts_bad, theta, x, x_int)),
      std::domain_error, expected_is_inf.str());
  ts_bad[0] = -inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts_bad, theta, x, x_int)),
      std::domain_error, "times");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts_bad, theta, x, x_int)),
      std::domain_error, expected_is_neg_inf.str());

  std::vector<double> theta_bad = theta;
  theta_bad[0] = inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts, theta_bad, x, x_int)),
      std::domain_error, "parameters and data");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts, theta_bad, x, x_int)),
      std::domain_error, expected_is_inf.str());
  theta_bad[0] = -inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts, theta_bad, x, x_int)),
      std::domain_error, "parameters and data");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0, t0, ts, theta_bad, x, x_int)),
      std::domain_error, expected_is_neg_inf.str());

  if (x.size() > 0) {
    std::vector<double> x_bad = x;
    x_bad[0] = inf;
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x_bad, x_int)),
        std::domain_error, "parameters and data");
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x_bad, x_int)),
        std::domain_error, expected_is_inf.str());
    x_bad[0] = -inf;
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x_bad, x_int)),
        std::domain_error, "parameters and data");
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(harm_osc, y0, t0, ts, theta, x_bad, x_int)),
        std::domain_error, expected_is_neg_inf.str());
  }
}

TEST(StanMathOde_integrate_ode_rk45, error_name) {
  using stan::math::integrate_ode_rk45;
  harm_osc_ode_data_fun harm_osc;

  std::vector<double> theta;
  theta.push_back(0.15);

  double t0 = 0;

  std::vector<double> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  std::vector<double> y0_bad;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(harm_osc, y0_bad, t0, ts, theta, x, x_int)),
      std::invalid_argument, "integrate_ode_rk45");
}
