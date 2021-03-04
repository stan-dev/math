#ifndef STAN_MATH_TEST_ODE_TEST_UTIL_HPP
#define STAN_MATH_TEST_ODE_TEST_UTIL_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <sstream>
#include <vector>

// calculates finite diffs for ode solver with varying parameters
template <typename F, typename solver_functor_t>
std::vector<Eigen::VectorXd> finite_diff_params(
    solver_functor_t& sol, const F& f, const double& t_in,
    const std::vector<double>& ts, const Eigen::VectorXd& y_in,
    const std::vector<double>& theta, const std::vector<double>& x,
    const std::vector<int>& x_int, const size_t& param_index,
    const double& diff) {
  std::stringstream msgs;
  std::vector<double> theta_ub(theta.size());
  std::vector<double> theta_lb(theta.size());
  for (size_t i = 0; i < theta.size(); i++) {
    if (i == param_index) {
      theta_ub[i] = theta[i] + diff;
      theta_lb[i] = theta[i] - diff;
    } else {
      theta_ub[i] = theta[i];
      theta_lb[i] = theta[i];
    }
  }

  std::vector<Eigen::VectorXd> ode_res_ub
      = sol(f, y_in, t_in, ts, &msgs, theta_ub, x, x_int);
  std::vector<Eigen::VectorXd> ode_res_lb
      = sol(f, y_in, t_in, ts, &msgs, theta_lb, x, x_int);

  std::vector<Eigen::VectorXd> results(ts.size());

  for (size_t i = 0; i < ts.size(); ++i) {
    results[i] = (ode_res_ub[i] - ode_res_lb[i]) / (2.0 * diff);
  }
  return results;
}

// calculates finite diffs for ode solver with varying initial positions
template <typename F, typename solver_functor_t>
std::vector<Eigen::VectorXd> finite_diff_initial_position(
    solver_functor_t& sol, const F& f, const double& t_in,
    const std::vector<double>& ts, const Eigen::VectorXd& y_in,
    const std::vector<double>& theta, const std::vector<double>& x,
    const std::vector<int>& x_int, const size_t& param_index,
    const double& diff) {
  std::stringstream msgs;
  Eigen::VectorXd y_in_ub(y_in);
  Eigen::VectorXd y_in_lb(y_in);
  for (size_t i = 0; i < y_in.size(); ++i) {
    if (i == param_index) {
      y_in_ub(i) = y_in(i) + diff;
      y_in_lb(i) = y_in(i) - diff;
    }
  }

  std::vector<Eigen::VectorXd> ode_res_ub
      = sol(f, y_in_ub, t_in, ts, &msgs, theta, x, x_int);
  std::vector<Eigen::VectorXd> ode_res_lb
      = sol(f, y_in_lb, t_in, ts, &msgs, theta, x, x_int);

  std::vector<Eigen::Matrix<double, -1, 1>> results(ts.size());

  for (size_t i = 0; i < ts.size(); ++i) {
    results[i] = (ode_res_ub[i] - ode_res_lb[i]) / (2 * diff);
  }
  return results;
}

// test ode solver with initial positions as doubles and parameters
// as vars against finite differences
template <typename F, typename solver_functor_t>
void test_ode_finite_diff_dv(solver_functor_t& sol, const F& f,
                             const double& t_in, const std::vector<double>& ts,
                             const Eigen::VectorXd& y_in,
                             const std::vector<double>& theta,
                             const std::vector<double>& x,
                             const std::vector<int>& x_int, const double& diff,
                             const double& diff2) {
  std::stringstream msgs;

  std::vector<std::vector<Eigen::VectorXd>> finite_diff_res(theta.size());
  for (size_t i = 0; i < theta.size(); i++) {
    finite_diff_res[i]
        = finite_diff_params(sol, f, t_in, ts, y_in, theta, x, x_int, i, diff);
  }

  std::vector<double> grads_eff;

  std::vector<stan::math::var> theta_v = stan::math::to_var(theta);

  std::vector<Eigen::Matrix<stan::math::var, -1, 1>> ode_res
      = sol(f, y_in, t_in, ts, &msgs, theta_v, x, x_int);

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = 0; j < y_in.size(); j++) {
      grads_eff.clear();
      ode_res[i][j].grad(theta_v, grads_eff);

      for (size_t k = 0; k < theta.size(); k++)
        EXPECT_NEAR(grads_eff[k], finite_diff_res[k][i][j], diff2)
            << "Gradient of ODE solver failed with initial positions"
            << " known and parameters unknown at time index " << i
            << ", equation index " << j << ", and parameter index: " << k;

      stan::math::set_zero_all_adjoints();
    }
  }
}

// test ode solver with initial positions as vars and parameters
// as doubles against finite differences
template <typename F, typename solver_functor_t>
void test_ode_finite_diff_vd(solver_functor_t& sol, const F& f,
                             const double& t_in, const std::vector<double>& ts,
                             const Eigen::VectorXd& y_in,
                             const std::vector<double>& theta,
                             const std::vector<double>& x,
                             const std::vector<int>& x_int, const double& diff,
                             const double& diff2) {
  std::stringstream msgs;

  std::vector<std::vector<Eigen::VectorXd>> finite_diff_res(y_in.size());
  for (size_t i = 0; i < y_in.size(); i++) {
    finite_diff_res[i] = finite_diff_initial_position(sol, f, t_in, ts, y_in,
                                                      theta, x, x_int, i, diff);
  }

  std::vector<double> grads_eff;

  Eigen::Matrix<stan::math::var, -1, 1> y_in_v = stan::math::to_var(y_in);

  std::vector<Eigen::Matrix<stan::math::var, -1, 1>> ode_res
      = sol(f, y_in_v, t_in, ts, &msgs, theta, x, x_int);

  std::vector<stan::math::var> y_vec(to_array_1d(y_in_v));

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = 0; j < y_in.size(); j++) {
      grads_eff.clear();
      ode_res[i][j].grad(y_vec, grads_eff);

      for (size_t k = 0; k < y_in.size(); k++)
        EXPECT_NEAR(grads_eff[k], finite_diff_res[k][i][j], diff2)
            << "Gradient of ode solver failed with initial positions"
            << " unknown and parameters known at time index " << i
            << ", equation index " << j << ", and parameter index: " << k;

      stan::math::set_zero_all_adjoints();
    }
  }
}

// test ode solver with initial positions as vars and parameters
// as vars against finite differences
template <typename F, typename solver_functor_t>
void test_ode_finite_diff_vv(solver_functor_t& sol, const F& f,
                             const double& t_in, const std::vector<double>& ts,
                             const Eigen::VectorXd& y_in,
                             const std::vector<double>& theta,
                             const std::vector<double>& x,
                             const std::vector<int>& x_int, const double& diff,
                             const double& diff2) {
  std::stringstream msgs;

  std::vector<std::vector<Eigen::VectorXd>> finite_diff_res_y(y_in.size());
  for (size_t i = 0; i < y_in.size(); i++) {
    finite_diff_res_y[i] = finite_diff_initial_position(
        sol, f, t_in, ts, y_in, theta, x, x_int, i, diff);
  }

  std::vector<std::vector<Eigen::VectorXd>> finite_diff_res_p(theta.size());
  for (size_t i = 0; i < theta.size(); i++)
    finite_diff_res_p[i]
        = finite_diff_params(sol, f, t_in, ts, y_in, theta, x, x_int, i, diff);

  std::vector<stan::math::var> vars(y_in.data(), y_in.data() + y_in.size());
  for (int i = 0; i < theta.size(); ++i) {
    vars.push_back(theta[i]);
  }
  Eigen::Matrix<stan::math::var, -1, 1> yv(y_in.size());
  for (int i = 0; i < yv.size(); ++i) {
    yv(i) = vars[i];
  }
  std::vector<stan::math::var> theta_v(vars.begin() + yv.size(), vars.end());

  std::vector<Eigen::Matrix<stan::math::var, -1, 1>> ode_res
      = sol(f, yv, t_in, ts, &msgs, theta_v, x, x_int);

  std::vector<double> grads_eff;
  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = 0; j < y_in.size(); j++) {
      grads_eff.clear();
      ode_res[i][j].grad(vars, grads_eff);

      for (size_t k = 0; k < theta.size(); k++)
        EXPECT_NEAR(grads_eff[k + y_in.size()], finite_diff_res_p[k][i][j],
                    diff2)
            << "Gradient of ode solver failed with initial positions"
            << " unknown and parameters unknown for param at time index " << i
            << ", equation index " << j << ", and parameter index: " << k;
      for (size_t k = 0; k < y_in.size(); k++)
        EXPECT_NEAR(grads_eff[k], finite_diff_res_y[k][i][j], diff2)
            << "Gradient of ode solver failed with initial positions"
            << " unknown and parameters known for initial position at time "
               "index "
            << i << ", equation index " << j << ", and parameter index: " << k;

      stan::math::set_zero_all_adjoints();
    }
  }
}

template <typename F, typename T1, typename T2, typename solver_functor_t>
void test_ode_error_conditions(solver_functor_t& sol, F& f, const double& t0,
                               const std::vector<double>& ts,
                               const Eigen::Matrix<T1, -1, 1>& y0,
                               const std::vector<T2>& theta,
                               const std::vector<double>& x,
                               const std::vector<int>& x_int) {
  std::stringstream msgs;

  ASSERT_NO_THROW((sol(f, y0, t0, ts, 0, theta, x, x_int)));
  ASSERT_EQ("", msgs.str());

  msgs.clear();
  Eigen::Matrix<T1, -1, 1> y0_bad;
  EXPECT_THROW_MSG((sol(f, y0_bad, t0, ts, &msgs, theta, x, x_int)),
                   std::invalid_argument, "initial state has size 0");
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  double t0_bad = ts[0] + 0.1;
  std::stringstream expected_msg;
  expected_msg << "initial time is " << t0_bad << ", but must be less than "
               << ts[0];
  EXPECT_THROW_MSG((sol(f, y0, t0_bad, ts, &msgs, theta, x, x_int)),
                   std::domain_error, expected_msg.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  std::vector<double> ts_bad;
  EXPECT_THROW_MSG((sol(f, y0, t0, ts_bad, &msgs, theta, x, x_int)),
                   std::invalid_argument, "times has size 0");
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  ts_bad.push_back(3);
  ts_bad.push_back(1);
  EXPECT_THROW_MSG((sol(f, y0, t0, ts_bad, &msgs, theta, x, x_int)),
                   std::domain_error, "times is not a valid sorted vector");
  EXPECT_EQ("", msgs.str());
}

template <typename F, typename T1, typename T2, typename solver_functor_t>
void test_ode_error_conditions_nan(solver_functor_t& sol, F& f,
                                   const double& t0,
                                   const std::vector<double>& ts,
                                   const Eigen::Matrix<T1, -1, 1>& y0,
                                   const std::vector<T2>& theta,
                                   const std::vector<double>& x,
                                   const std::vector<int>& x_int) {
  std::stringstream msgs;
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::stringstream expected_is_nan;
  expected_is_nan << "is " << nan;

  ASSERT_NO_THROW((sol(f, y0, t0, ts, 0, theta, x, x_int)));
  ASSERT_EQ("", msgs.str());

  msgs.clear();
  Eigen::Matrix<T1, -1, 1> y0_bad = y0;
  y0_bad[0] = nan;
  EXPECT_THROW_MSG((sol(f, y0_bad, t0, ts, &msgs, theta, x, x_int)),
                   std::domain_error, "initial state");
  EXPECT_THROW_MSG((sol(f, y0_bad, t0, ts, &msgs, theta, x, x_int)),
                   std::domain_error, expected_is_nan.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  double t0_bad = nan;
  EXPECT_THROW_MSG((sol(f, y0, t0_bad, ts, &msgs, theta, x, x_int)),
                   std::domain_error, "initial time");
  EXPECT_THROW_MSG((sol(f, y0, t0_bad, ts, &msgs, theta, x, x_int)),
                   std::domain_error, expected_is_nan.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  std::vector<double> ts_bad = ts;
  ts_bad[0] = nan;
  EXPECT_THROW_MSG((sol(f, y0, t0, ts_bad, &msgs, theta, x, x_int)),
                   std::domain_error, "times");
  EXPECT_THROW_MSG((sol(f, y0, t0, ts_bad, &msgs, theta, x, x_int)),
                   std::domain_error, expected_is_nan.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  std::vector<T2> theta_bad = theta;
  theta_bad[0] = nan;
  EXPECT_THROW_MSG((sol(f, y0, t0, ts, &msgs, theta_bad, x, x_int)),
                   std::domain_error, "ode parameters and data");
  EXPECT_THROW_MSG((sol(f, y0, t0, ts, &msgs, theta_bad, x, x_int)),
                   std::domain_error, expected_is_nan.str());
  EXPECT_EQ("", msgs.str());

  if (x.size() > 0) {
    msgs.clear();
    std::vector<double> x_bad = x;
    x_bad[0] = nan;
    EXPECT_THROW_MSG((sol(f, y0, t0, ts, &msgs, theta, x_bad, x_int)),
                     std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG((sol(f, y0, t0, ts, &msgs, theta, x_bad, x_int)),
                     std::domain_error, expected_is_nan.str());
    EXPECT_EQ("", msgs.str());
  }
}

template <typename F, typename T1, typename T2, typename solver_functor_t>
void test_ode_error_conditions_inf(solver_functor_t& sol, F& f,
                                   const double& t0,
                                   const std::vector<double>& ts,
                                   const Eigen::Matrix<T1, -1, 1>& y0,
                                   const std::vector<T2>& theta,
                                   const std::vector<double>& x,
                                   const std::vector<int>& x_int) {
  std::stringstream msgs;
  double inf = std::numeric_limits<double>::infinity();
  std::stringstream expected_is_inf;
  expected_is_inf << "is " << inf;
  std::stringstream expected_is_neg_inf;
  expected_is_neg_inf << "is " << -inf;

  ASSERT_NO_THROW((sol(f, y0, t0, ts, 0, theta, x, x_int)));
  ASSERT_EQ("", msgs.str());

  msgs.clear();
  Eigen::Matrix<T1, -1, 1> y0_bad = y0;
  y0_bad[0] = inf;
  EXPECT_THROW_MSG((sol(f, y0_bad, t0, ts, &msgs, theta, x, x_int)),
                   std::domain_error, "initial state");
  EXPECT_THROW_MSG((sol(f, y0_bad, t0, ts, &msgs, theta, x, x_int)),
                   std::domain_error, expected_is_inf.str());
  y0_bad[0] = -inf;
  EXPECT_THROW_MSG((sol(f, y0_bad, t0, ts, &msgs, theta, x, x_int)),
                   std::domain_error, "initial state");
  EXPECT_THROW_MSG((sol(f, y0_bad, t0, ts, &msgs, theta, x, x_int)),
                   std::domain_error, expected_is_neg_inf.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  double t0_bad = inf;
  EXPECT_THROW_MSG((sol(f, y0, t0_bad, ts, &msgs, theta, x, x_int)),
                   std::domain_error, "initial time");
  EXPECT_THROW_MSG((sol(f, y0, t0_bad, ts, &msgs, theta, x, x_int)),
                   std::domain_error, expected_is_inf.str());
  t0_bad = -inf;
  EXPECT_THROW_MSG((sol(f, y0, t0_bad, ts, &msgs, theta, x, x_int)),
                   std::domain_error, "initial time");
  EXPECT_THROW_MSG((sol(f, y0, t0_bad, ts, &msgs, theta, x, x_int)),
                   std::domain_error, expected_is_neg_inf.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  std::vector<double> ts_bad = ts;
  ts_bad[0] = inf;
  EXPECT_THROW_MSG((sol(f, y0, t0, ts_bad, &msgs, theta, x, x_int)),
                   std::domain_error, "times");
  EXPECT_THROW_MSG((sol(f, y0, t0, ts_bad, &msgs, theta, x, x_int)),
                   std::domain_error, expected_is_inf.str());
  ts_bad[0] = -inf;
  EXPECT_THROW_MSG((sol(f, y0, t0, ts_bad, &msgs, theta, x, x_int)),
                   std::domain_error, "times");
  EXPECT_THROW_MSG((sol(f, y0, t0, ts_bad, &msgs, theta, x, x_int)),
                   std::domain_error, expected_is_neg_inf.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  std::vector<T2> theta_bad = theta;
  theta_bad[0] = inf;
  EXPECT_THROW_MSG((sol(f, y0, t0, ts, &msgs, theta_bad, x, x_int)),
                   std::domain_error, "ode parameters and data");
  EXPECT_THROW_MSG((sol(f, y0, t0, ts, &msgs, theta_bad, x, x_int)),
                   std::domain_error, expected_is_inf.str());
  theta_bad[0] = -inf;
  EXPECT_THROW_MSG((sol(f, y0, t0, ts, &msgs, theta_bad, x, x_int)),
                   std::domain_error, "ode parameters and data");
  EXPECT_THROW_MSG((sol(f, y0, t0, ts, &msgs, theta_bad, x, x_int)),
                   std::domain_error, expected_is_neg_inf.str());
  EXPECT_EQ("", msgs.str());

  if (x.size() > 0) {
    msgs.clear();
    std::vector<double> x_bad = x;
    x_bad[0] = inf;
    EXPECT_THROW_MSG((sol(f, y0, t0, ts, &msgs, theta, x_bad, x_int)),
                     std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG((sol(f, y0, t0, ts, &msgs, theta, x_bad, x_int)),
                     std::domain_error, expected_is_inf.str());
    x_bad[0] = -inf;
    EXPECT_THROW_MSG((sol(f, y0, t0, ts, &msgs, theta, x_bad, x_int)),
                     std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG((sol(f, y0, t0, ts, &msgs, theta, x_bad, x_int)),
                     std::domain_error, expected_is_neg_inf.str());
    EXPECT_EQ("", msgs.str());
  }
}

template <typename F, typename solver_functor_t>
void test_ode(solver_functor_t& sol, const F& f, const double& t_in,
              const std::vector<double>& ts, const Eigen::VectorXd& y_in,
              const std::vector<double>& theta, const std::vector<double>& x,
              const std::vector<int>& x_int, const double& diff,
              const double& diff2) {
  using stan::math::to_var;

  test_ode_finite_diff_vd(sol, f, t_in, ts, y_in, theta, x, x_int, diff, diff2);
  test_ode_finite_diff_dv(sol, f, t_in, ts, y_in, theta, x, x_int, diff, diff2);
  test_ode_finite_diff_vv(sol, f, t_in, ts, y_in, theta, x, x_int, diff, diff2);

  // vd
  {
    Eigen::Matrix<stan::math::var, -1, 1> y_var = to_var(y_in);
    test_ode_error_conditions(sol, f, t_in, ts, y_var, theta, x, x_int);
    test_ode_error_conditions_nan(sol, f, t_in, ts, y_var, theta, x, x_int);
  }

  // dv
  {
    std::vector<stan::math::var> theta_var = to_var(theta);
    test_ode_error_conditions(sol, f, t_in, ts, y_in, theta_var, x, x_int);
    test_ode_error_conditions_nan(sol, f, t_in, ts, y_in, theta_var, x, x_int);
  }

  // vv
  {
    Eigen::Matrix<stan::math::var, -1, 1> y_var = to_var(y_in);
    std::vector<stan::math::var> theta_var = to_var(theta);
    test_ode_error_conditions(sol, f, t_in, ts, y_var, theta_var, x, x_int);
    test_ode_error_conditions_nan(sol, f, t_in, ts, y_var, theta_var, x, x_int);
    test_ode_error_conditions_inf(sol, f, t_in, ts, y_var, theta_var, x, x_int);
  }
}

#endif
