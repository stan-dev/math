#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/rev/functor/coupled_mm.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

template <typename F>
void sho_value_test(F harm_osc, std::vector<double>& y0, double t0,
                    std::vector<double>& ts, std::vector<double>& theta,
                    std::vector<double>& x, std::vector<int>& x_int) {
  std::vector<std::vector<double> > ode_res_vd
      = stan::math::integrate_ode_adams(harm_osc, y0, t0, ts, theta, x, x_int,
                                        0, 1e-8, 1e-10, 1e6);

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

TEST(StanMathOde_integrate_ode_adams, harmonic_oscillator) {
  sho_test(0.0);
  sho_test(1.0);
  sho_test(-1.0);

  sho_data_test(0.0);
  sho_data_test(1.0);
  sho_data_test(-1.0);
}

TEST(StanMathOde_integrate_ode_adams, error_conditions) {
  using stan::math::integrate_ode_adams;
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

  ASSERT_NO_THROW(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x, x_int, 0,
                                      1e-8, 1e-10, 1e6));

  std::vector<double> y0_bad;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0_bad, t0, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::invalid_argument, "initial state has size 0");

  double t0_bad = 2.0;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0_bad, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error,
                   "initial time is 2, but must be less than 0.1");

  std::vector<double> ts_bad;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts_bad, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::invalid_argument, "times has size 0");

  ts_bad.push_back(3);
  ts_bad.push_back(1);
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts_bad, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "times is not a valid sorted vector");

  // TODO(carpenter): g++6 failure
  std::vector<double> theta_bad;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta_bad, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::out_of_range, "vector");

  // TODO(carpenter): g++6 failure
  std::vector<double> x_bad;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x_bad,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::out_of_range, "vector");

  // TODO(carpenter): g++6 failure
  std::vector<int> x_int_bad;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x,
                                       x_int_bad, 0, 1e-8, 1e-10, 1e6),
                   std::out_of_range, "vector");

  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x, x_int, 0,
                                       -1, 1e-6, 10),
                   std::domain_error, "relative_tolerance");

  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x, x_int, 0,
                                       1e-6, -1, 10),
                   std::domain_error, "absolute_tolerance");

  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x, x_int, 0,
                                       1e-6, 1e-6, -1),
                   std::domain_error, "max_num_steps");
}

TEST(StanMathOde_integrate_ode_adams, error_conditions_nan) {
  using stan::math::integrate_ode_adams;
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

  ASSERT_NO_THROW(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x, x_int, 0,
                                      1e-8, 1e-10, 1e6));

  double nan = std::numeric_limits<double>::quiet_NaN();
  std::stringstream expected_is_nan;
  expected_is_nan << "is " << nan;

  std::vector<double> y0_bad = y0;
  y0_bad[0] = nan;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0_bad, t0, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "initial state");
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0_bad, t0, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, expected_is_nan.str());

  double t0_bad = nan;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0_bad, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "initial time");
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0_bad, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, expected_is_nan.str());

  std::vector<double> ts_bad = ts;
  ts_bad[0] = nan;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts_bad, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "times");
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts_bad, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, expected_is_nan.str());

  std::vector<double> theta_bad = theta;
  theta_bad[0] = nan;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta_bad, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "ode parameters and data");
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta_bad, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, expected_is_nan.str());

  if (x.size() > 0) {
    std::vector<double> x_bad = x;
    x_bad[0] = nan;
    EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x_bad,
                                         x_int, 0, 1e-8, 1e-10, 1e6),
                     std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x_bad,
                                         x_int, 0, 1e-8, 1e-10, 1e6),
                     std::domain_error, expected_is_nan.str());
  }
}

TEST(StanMathOde_integrate_ode_adams, error_conditions_inf) {
  using stan::math::integrate_ode_adams;
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

  ASSERT_NO_THROW(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x, x_int, 0,
                                      1e-8, 1e-10, 1e6));

  double inf = std::numeric_limits<double>::infinity();
  std::vector<double> y0_bad = y0;
  y0_bad[0] = inf;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0_bad, t0, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "initial state");
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0_bad, t0, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, expected_is_inf.str());

  y0_bad[0] = -inf;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0_bad, t0, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "initial state");
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0_bad, t0, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, expected_is_neg_inf.str());

  double t0_bad = inf;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0_bad, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "initial time");
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0_bad, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, expected_is_inf.str());
  t0_bad = -inf;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0_bad, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "initial time");
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0_bad, ts, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, expected_is_neg_inf.str());

  std::vector<double> ts_bad = ts;
  ts_bad.push_back(inf);
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts_bad, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "times");
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts_bad, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, expected_is_inf.str());
  ts_bad[0] = -inf;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts_bad, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "times");
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts_bad, theta, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, expected_is_neg_inf.str());

  std::vector<double> theta_bad = theta;
  theta_bad[0] = inf;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta_bad, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "ode parameters and data");
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta_bad, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, expected_is_inf.str());
  theta_bad[0] = -inf;
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta_bad, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, "ode parameters and data");
  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta_bad, x,
                                       x_int, 0, 1e-8, 1e-10, 1e6),
                   std::domain_error, expected_is_neg_inf.str());

  if (x.size() > 0) {
    std::vector<double> x_bad = x;
    x_bad[0] = inf;
    EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x_bad,
                                         x_int, 0, 1e-8, 1e-10, 1e6),
                     std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x_bad,
                                         x_int, 0, 1e-8, 1e-10, 1e6),
                     std::domain_error, expected_is_inf.str());
    x_bad[0] = -inf;
    EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x_bad,
                                         x_int, 0, 1e-8, 1e-10, 1e6),
                     std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x_bad,
                                         x_int, 0, 1e-8, 1e-10, 1e6),
                     std::domain_error, expected_is_neg_inf.str());
  }
}

TEST(StanMathOde_integrate_ode_adams, error_conditions_bad_ode) {
  using stan::math::integrate_ode_adams;
  harm_osc_ode_wrong_size_1_fun harm_osc;

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

  std::string error_msg
      = "cvodes_integrator: dy_dt (3) and states (2) must match in size";

  EXPECT_THROW_MSG(integrate_ode_adams(harm_osc, y0, t0, ts, theta, x, x_int, 0,
                                       1e-8, 1e-10, 1e6),
                   std::invalid_argument, error_msg);
}

TEST(StanMathOde_integrate_ode_adams, too_much_work) {
  coupled_mm_ode_fun f_;

  // initial value and parameters from model definition
  std::vector<double> y0(2);
  y0[0] = 1E5;
  y0[1] = 1E-1;

  double t0 = 0;

  std::vector<double> ts_long;
  ts_long.push_back(1E10);

  std::vector<double> ts_short;
  ts_short.push_back(1);

  std::vector<double> theta(4);

  theta[0] = 1.0;
  theta[1] = 0.5;
  theta[2] = 0.5;
  theta[3] = 0.1;

  std::vector<double> data;

  std::vector<int> data_int;

  EXPECT_THROW_MSG(
      stan::math::integrate_ode_adams(f_, y0, t0, ts_long, theta, data,
                                      data_int, 0, 1E-6, 1E-6, 100),
      std::domain_error,
      "integrate_ode_adams:  Failed to integrate to next output time");

  EXPECT_NO_THROW(stan::math::integrate_ode_bdf(
      f_, y0, t0, ts_short, theta, data, data_int, 0, 1E-6, 1E-6, 100));
}
