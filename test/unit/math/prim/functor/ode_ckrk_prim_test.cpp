#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <stan/math/prim/functor/ode_ckrk.hpp>
#include <boost/numeric/odeint.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/functor/lorenz.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>

using stan::math::ode_ckrk;
using stan::math::ode_ckrk_tol;
using stan::math::to_vector;

TEST(StanMathOde_ode_ckrk, harmonic_oscillator_value) {
  using stan::math::ode_ckrk;
  std::vector<double> theta{0.15};
  Eigen::VectorXd y0(2);
  y0 << 1.0, 0.0;
  std::vector<double> ts(100);
  std::vector<double> x_r;
  std::vector<int> x_i;
  double t0;
  harm_osc_ode_fun_eigen harm_osc;

  {
    t0 = 0.0;
    for (int i = 0; i < ts.size(); ++i) {
      ts[i] = t0 + 0.1 * (i + 1);
    }
    auto res = ode_ckrk(harm_osc, y0, t0, ts, 0, theta, x_r, x_i);
    EXPECT_NEAR(0.995029, res[0][0], 1e-5);
    EXPECT_NEAR(-0.0990884, res[0][1], 1e-5);

    EXPECT_NEAR(-0.421907, res[99][0], 1e-5);
    EXPECT_NEAR(0.246407, res[99][1], 1e-5);
  }

  {
    t0 = -1.0;
    for (int i = 0; i < ts.size(); ++i) {
      ts[i] = t0 + 0.1 * (i + 1);
    }
    auto res = ode_ckrk(harm_osc, y0, t0, ts, 0, theta, x_r, x_i);
    EXPECT_NEAR(0.995029, res[0][0], 1e-5);
    EXPECT_NEAR(-0.0990884, res[0][1], 1e-5);

    EXPECT_NEAR(-0.421907, res[99][0], 1e-5);
    EXPECT_NEAR(0.246407, res[99][1], 1e-5);
  }

  {
    t0 = 1.0;
    for (int i = 0; i < ts.size(); ++i) {
      ts[i] = t0 + 0.1 * (i + 1);
    }
    auto res = ode_ckrk(harm_osc, y0, t0, ts, 0, theta, x_r, x_i);
    EXPECT_NEAR(0.995029, res[0][0], 1e-5);
    EXPECT_NEAR(-0.0990884, res[0][1], 1e-5);

    EXPECT_NEAR(-0.421907, res[99][0], 1e-5);
    EXPECT_NEAR(0.246407, res[99][1], 1e-5);    
  }
}

TEST(StanMathOde_ode_ckrk, harmonic_oscillator_data_value) {
  using stan::math::ode_ckrk;
  std::vector<double> theta{0.15};
  Eigen::VectorXd y0(2);
  y0 << 1.0, 0.0;
  std::vector<double> ts(100);
  std::vector<double> x_r(3, 1);
  std::vector<int> x_i(2, 0);
  double t0;
  harm_osc_ode_data_fun_eigen harm_osc;

  {
    t0 = 0.0;
    for (int i = 0; i < ts.size(); ++i) {
      ts[i] = t0 + 0.1 * (i + 1);
    }
    auto res = ode_ckrk(harm_osc, y0, t0, ts, 0, theta, x_r, x_i);
    EXPECT_NEAR(0.995029, res[0][0], 1e-5);
    EXPECT_NEAR(-0.0990884, res[0][1], 1e-5);

    EXPECT_NEAR(-0.421907, res[99][0], 1e-5);
    EXPECT_NEAR(0.246407, res[99][1], 1e-5);
  }

  {
    t0 = -1.0;
    for (int i = 0; i < ts.size(); ++i) {
      ts[i] = t0 + 0.1 * (i + 1);
    }
    auto res = ode_ckrk(harm_osc, y0, t0, ts, 0, theta, x_r, x_i);
    EXPECT_NEAR(0.995029, res[0][0], 1e-5);
    EXPECT_NEAR(-0.0990884, res[0][1], 1e-5);

    EXPECT_NEAR(-0.421907, res[99][0], 1e-5);
    EXPECT_NEAR(0.246407, res[99][1], 1e-5);
  }

  {
    t0 = 1.0;
    for (int i = 0; i < ts.size(); ++i) {
      ts[i] = t0 + 0.1 * (i + 1);
    }
    auto res = ode_ckrk(harm_osc, y0, t0, ts, 0, theta, x_r, x_i);
    EXPECT_NEAR(0.995029, res[0][0], 1e-5);
    EXPECT_NEAR(-0.0990884, res[0][1], 1e-5);

    EXPECT_NEAR(-0.421907, res[99][0], 1e-5);
    EXPECT_NEAR(0.246407, res[99][1], 1e-5);    
  }
}

TEST(StanMathOde_ode_ckrk, error_conditions) {
  std::ostream* msgs = nullptr;
  harm_osc_ode_data_fun_eigen harm_osc;

  std::vector<double> theta;
  theta.push_back(0.15);

  Eigen::VectorXd y0(2);
  y0 << 1.0, 0.0;

  double t0 = 0;

  std::vector<double> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  ASSERT_NO_THROW((ode_ckrk(harm_osc, y0, t0, ts, msgs, theta, x, x_int)));

  Eigen::VectorXd y0_bad;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0_bad, t0, ts, msgs, theta, x, x_int)),
                   std::invalid_argument, "initial state has size 0");

  double t0_bad = 2.0;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0_bad, ts, msgs, theta, x, x_int)),
                   std::domain_error,
                   "initial time is 2, but must be less than 0.1");

  std::vector<double> ts_bad;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts_bad, msgs, theta, x, x_int)),
                   std::invalid_argument, "times has size 0");

  ts_bad.push_back(3);
  ts_bad.push_back(1);
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts_bad, msgs, theta, x, x_int)),
                   std::domain_error, "times is not a valid sorted vector");

  std::vector<double> theta_bad;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts, msgs, theta_bad, x, x_int)),
                   std::out_of_range, "vector");

  std::vector<double> x_bad;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts, msgs, theta, x_bad, x_int)),
                   std::out_of_range, "vector");

  std::vector<int> x_int_bad;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts, msgs, theta, x, x_int_bad)),
                   std::out_of_range, "vector");

  EXPECT_THROW_MSG(
      (ode_ckrk_tol(harm_osc, y0, t0, ts, -1, 1e-6, 10, msgs, theta, x, x_int)),
      std::domain_error, "relative_tolerance");

  EXPECT_THROW_MSG(
      (ode_ckrk_tol(harm_osc, y0, t0, ts, 1e-6, -1, 10, msgs, theta, x, x_int)),
      std::domain_error, "absolute_tolerance");

  EXPECT_THROW_MSG((ode_ckrk_tol(harm_osc, y0, t0, ts, 1e-6, 1e-6, -1, msgs,
                                 theta, x, x_int)),
                   std::domain_error, "max_num_steps");
}

TEST(StanMathOde_ode_ckrk, error_conditions_nan) {
  std::ostream* msgs = nullptr;
  harm_osc_ode_data_fun_eigen harm_osc;

  std::vector<double> theta;
  theta.push_back(0.15);

  Eigen::VectorXd y0(2);
  y0 << 1.0, 2.0;

  double t0 = 0;

  std::vector<double> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  ASSERT_NO_THROW((ode_ckrk(harm_osc, y0, t0, ts, msgs, theta, x, x_int)));

  double nan = std::numeric_limits<double>::quiet_NaN();
  std::stringstream expected_is_nan;
  expected_is_nan << "is " << nan;

  Eigen::VectorXd y0_bad = y0;
  y0_bad[0] = nan;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0_bad, t0, ts, msgs, theta, x, x_int)),
                   std::domain_error, "initial state");
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0_bad, t0, ts, msgs, theta, x, x_int)),
                   std::domain_error, expected_is_nan.str());

  double t0_bad = nan;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0_bad, ts, msgs, theta, x, x_int)),
                   std::domain_error, "initial time");
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0_bad, ts, msgs, theta, x, x_int)),
                   std::domain_error, expected_is_nan.str());

  std::vector<double> ts_bad = ts;
  ts_bad[0] = nan;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts_bad, msgs, theta, x, x_int)),
                   std::domain_error, "times");
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts_bad, msgs, theta, x, x_int)),
                   std::domain_error, expected_is_nan.str());

  std::vector<double> theta_bad = theta;
  theta_bad[0] = nan;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts, msgs, theta_bad, x, x_int)),
                   std::domain_error, "parameters and data");
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts, msgs, theta_bad, x, x_int)),
                   std::domain_error, expected_is_nan.str());

  if (x.size() > 0) {
    std::vector<double> x_bad = x;
    x_bad[0] = nan;
    EXPECT_THROW_MSG(
        (ode_ckrk(harm_osc, y0, t0, ts, msgs, theta, x_bad, x_int)),
        std::domain_error, "parameters and data");
    EXPECT_THROW_MSG(
        (ode_ckrk(harm_osc, y0, t0, ts, msgs, theta, x_bad, x_int)),
        std::domain_error, expected_is_nan.str());
  }
}

TEST(StanMathOde_ode_ckrk, error_conditions_inf) {
  std::ostream* msgs = nullptr;
  harm_osc_ode_data_fun_eigen harm_osc;
  std::stringstream expected_is_inf;
  expected_is_inf << "is " << std::numeric_limits<double>::infinity();
  std::stringstream expected_is_neg_inf;
  expected_is_neg_inf << "is " << -std::numeric_limits<double>::infinity();

  std::vector<double> theta;
  theta.push_back(0.15);

  Eigen::VectorXd y0(2);
  y0 << 1.0, 2.0;

  double t0 = 0;

  std::vector<double> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  ASSERT_NO_THROW((ode_ckrk(harm_osc, y0, t0, ts, msgs, theta, x, x_int)));

  double inf = std::numeric_limits<double>::infinity();
  Eigen::VectorXd y0_bad = y0;
  y0_bad[0] = inf;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0_bad, t0, ts, msgs, theta, x, x_int)),
                   std::domain_error, "initial state");
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0_bad, t0, ts, msgs, theta, x, x_int)),
                   std::domain_error, expected_is_inf.str());

  y0_bad[0] = -inf;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0_bad, t0, ts, msgs, theta, x, x_int)),
                   std::domain_error, "initial state");
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0_bad, t0, ts, msgs, theta, x, x_int)),
                   std::domain_error, expected_is_neg_inf.str());

  double t0_bad = inf;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0_bad, ts, msgs, theta, x, x_int)),
                   std::domain_error, "initial time");
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0_bad, ts, msgs, theta, x, x_int)),
                   std::domain_error, expected_is_inf.str());
  t0_bad = -inf;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0_bad, ts, msgs, theta, x, x_int)),
                   std::domain_error, "initial time");
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0_bad, ts, msgs, theta, x, x_int)),
                   std::domain_error, expected_is_neg_inf.str());

  std::vector<double> ts_bad = ts;
  ts_bad.push_back(inf);
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts_bad, msgs, theta, x, x_int)),
                   std::domain_error, "times");
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts_bad, msgs, theta, x, x_int)),
                   std::domain_error, expected_is_inf.str());
  ts_bad[0] = -inf;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts_bad, msgs, theta, x, x_int)),
                   std::domain_error, "times");
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts_bad, msgs, theta, x, x_int)),
                   std::domain_error, expected_is_neg_inf.str());

  std::vector<double> theta_bad = theta;
  theta_bad[0] = inf;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts, msgs, theta_bad, x, x_int)),
                   std::domain_error, "parameters and data");
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts, msgs, theta_bad, x, x_int)),
                   std::domain_error, expected_is_inf.str());
  theta_bad[0] = -inf;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts, msgs, theta_bad, x, x_int)),
                   std::domain_error, "parameters and data");
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0, t0, ts, msgs, theta_bad, x, x_int)),
                   std::domain_error, expected_is_neg_inf.str());

  if (x.size() > 0) {
    std::vector<double> x_bad = x;
    x_bad[0] = inf;
    EXPECT_THROW_MSG(
        (ode_ckrk(harm_osc, y0, t0, ts, msgs, theta, x_bad, x_int)),
        std::domain_error, "parameters and data");
    EXPECT_THROW_MSG(
        (ode_ckrk(harm_osc, y0, t0, ts, msgs, theta, x_bad, x_int)),
        std::domain_error, expected_is_inf.str());
    x_bad[0] = -inf;
    EXPECT_THROW_MSG(
        (ode_ckrk(harm_osc, y0, t0, ts, msgs, theta, x_bad, x_int)),
        std::domain_error, "parameters and data");
    EXPECT_THROW_MSG(
        (ode_ckrk(harm_osc, y0, t0, ts, msgs, theta, x_bad, x_int)),
        std::domain_error, expected_is_neg_inf.str());
  }
}

TEST(StanMathOde_ode_ckrk, error_name) {
  std::ostream* msgs = nullptr;
  harm_osc_ode_data_fun_eigen harm_osc;

  std::vector<double> theta;
  theta.push_back(0.15);

  double t0 = 0;

  std::vector<double> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  Eigen::VectorXd y0_bad;
  EXPECT_THROW_MSG((ode_ckrk(harm_osc, y0_bad, t0, ts, msgs, theta, x, x_int)),
                   std::invalid_argument, "ode_ckrk");
}
