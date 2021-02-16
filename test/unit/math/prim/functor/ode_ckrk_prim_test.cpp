#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <stan/math/prim/functor/ode_ckrk.hpp>
#include <boost/numeric/odeint.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/functor/lorenz.hpp>
#include <test/unit/math/prim/functor/ode_test_functors.hpp>
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

// errors
TEST(ode_ckrk_prim, y0_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  Eigen::VectorXd y0inf(1);
  Eigen::VectorXd y0NaN(1);
  Eigen::VectorXd y0_empty;
  y0NaN << stan::math::NOT_A_NUMBER;
  y0inf << stan::math::INFTY;
  int t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, a));

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0inf, t0, ts, nullptr, a),
      std::domain_error);

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0NaN, t0, ts, nullptr, a),
      std::domain_error);

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0_empty, t0, ts, nullptr, a),
      std::invalid_argument);
}

TEST(ode_ckrk_prim, t0_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  double t0inf = stan::math::INFTY;
  double t0NaN = stan::math::NOT_A_NUMBER;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, a));

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0inf, ts, nullptr, a),
      std::domain_error);

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0NaN, ts, nullptr, a),
      std::domain_error);
}

TEST(ode_ckrk_prim, ts_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};
  std::vector<double> ts_repeat = {0.45, 0.45};
  std::vector<double> ts_lots = {0.45, 0.45, 1.1, 1.1, 2.0};
  std::vector<double> ts_empty = {};
  std::vector<double> ts_early = {-0.45, 0.2};
  std::vector<double> ts_decreasing = {0.45, 0.2};
  std::vector<double> tsinf = {stan::math::INFTY, 1.1};
  std::vector<double> tsNaN = {0.45, stan::math::NOT_A_NUMBER};

  double a = 1.5;

  std::vector<Eigen::VectorXd> out;
  EXPECT_NO_THROW(out = stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts,
                                             nullptr, a));
  EXPECT_EQ(out.size(), ts.size());

  EXPECT_NO_THROW(out = stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0,
                                             ts_repeat, nullptr, a));
  EXPECT_EQ(out.size(), ts_repeat.size());
  EXPECT_MATRIX_FLOAT_EQ(out[0], out[1]);

  EXPECT_NO_THROW(out = stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0,
                                             ts_lots, nullptr, a));
  EXPECT_EQ(out.size(), ts_lots.size());

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts_empty, nullptr, a),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts_early, nullptr, a),
      std::domain_error);

  EXPECT_THROW(stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0,
                                    ts_decreasing, nullptr, a),
               std::domain_error);

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, tsinf, nullptr, a),
      std::domain_error);

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, tsNaN, nullptr, a),
      std::domain_error);
}

TEST(ode_ckrk_prim, one_arg_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;
  double ainf = stan::math::INFTY;
  double aNaN = stan::math::NOT_A_NUMBER;

  std::vector<double> va = {a};
  std::vector<double> vainf = {ainf};
  std::vector<double> vaNaN = {aNaN};

  Eigen::VectorXd ea(1);
  ea << a;
  Eigen::VectorXd eainf(1);
  eainf << ainf;
  Eigen::VectorXd eaNaN(1);
  eaNaN << aNaN;

  std::vector<std::vector<double>> vva = {va};
  std::vector<std::vector<double>> vvainf = {vainf};
  std::vector<std::vector<double>> vvaNaN = {vaNaN};

  std::vector<Eigen::VectorXd> vea = {ea};
  std::vector<Eigen::VectorXd> veainf = {eainf};
  std::vector<Eigen::VectorXd> veaNaN = {eaNaN};

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, a));

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, ainf),
      std::domain_error);

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, aNaN),
      std::domain_error);

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, va));

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, vainf),
      std::domain_error);

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, vaNaN),
      std::domain_error);

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, ea));

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, eainf),
      std::domain_error);

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, eaNaN),
      std::domain_error);

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, vva));

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, vvainf),
      std::domain_error);

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, vvaNaN),
      std::domain_error);

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, vea));

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, veainf),
      std::domain_error);

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, veaNaN),
      std::domain_error);
}

TEST(ode_ckrk_prim, two_arg_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;
  double ainf = stan::math::INFTY;
  double aNaN = stan::math::NOT_A_NUMBER;

  std::vector<double> va = {a};
  std::vector<double> vainf = {ainf};
  std::vector<double> vaNaN = {aNaN};

  Eigen::VectorXd ea(1);
  ea << a;
  Eigen::VectorXd eainf(1);
  eainf << ainf;
  Eigen::VectorXd eaNaN(1);
  eaNaN << aNaN;

  std::vector<std::vector<double>> vva = {va};
  std::vector<std::vector<double>> vvainf = {vainf};
  std::vector<std::vector<double>> vvaNaN = {vaNaN};

  std::vector<Eigen::VectorXd> vea = {ea};
  std::vector<Eigen::VectorXd> veainf = {eainf};
  std::vector<Eigen::VectorXd> veaNaN = {eaNaN};

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr, a, a));

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr, a, ainf),
      std::domain_error);

  EXPECT_THROW(
      stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr, a, aNaN),
      std::domain_error);

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr, a, va));

  EXPECT_THROW(stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr,
                                    a, vainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr,
                                    a, vaNaN),
               std::domain_error);

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr, a, ea));

  EXPECT_THROW(stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr,
                                    a, eainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr,
                                    a, eaNaN),
               std::domain_error);

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr, a, vva));

  EXPECT_THROW(stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr,
                                    a, vvainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr,
                                    a, vvaNaN),
               std::domain_error);

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr, a, vea));

  EXPECT_THROW(stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr,
                                    a, veainf),
               std::domain_error);

  EXPECT_THROW(stan::math::ode_ckrk(stan::test::Cos2Arg(), y0, t0, ts, nullptr,
                                    a, veaNaN),
               std::domain_error);
}

TEST(ode_ckrk_prim, rhs_wrong_size_errors) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;

  EXPECT_NO_THROW(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, a));

  EXPECT_THROW(stan::math::ode_ckrk(stan::test::CosArgWrongSize(), y0, t0, ts,
                                    nullptr, a),
               std::invalid_argument);
}

TEST(ode_ckrk_prim, error_name) {
  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double ainf = stan::math::INFTY;

  EXPECT_THROW_MSG(
      stan::math::ode_ckrk(stan::test::CosArg1(), y0, t0, ts, nullptr, ainf),
      std::domain_error, "ode_ckrk");
}
