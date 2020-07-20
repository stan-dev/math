#include <stan/math/rev/core.hpp>
#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <sstream>
#include <vector>

// calculates finite diffs for integrate_ode_rk45 with varying parameters
template <typename F>
std::vector<std::vector<double>> finite_diff_params(
    const F& f, const double& t_in, const std::vector<double>& ts,
    const std::vector<double>& y_in, const std::vector<double>& theta,
    const std::vector<double>& x, const std::vector<int>& x_int,
    const size_t& param_index, const double& diff) {
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

  std::vector<std::vector<double>> ode_res_ub;
  std::vector<std::vector<double>> ode_res_lb;

  ode_res_ub = stan::math::integrate_ode_rk45(f, y_in, t_in, ts, theta_ub, x,
                                              x_int, &msgs);
  ode_res_lb = stan::math::integrate_ode_rk45(f, y_in, t_in, ts, theta_lb, x,
                                              x_int, &msgs);

  std::vector<std::vector<double>> results(ts.size());

  for (size_t i = 0; i < ode_res_ub.size(); i++)
    for (size_t j = 0; j < ode_res_ub[j].size(); j++)
      results[i].push_back((ode_res_ub[i][j] - ode_res_lb[i][j]) / (2 * diff));
  return results;
}

// calculates finite diffs for integrate_ode_rk45 with varying initial positions
template <typename F>
std::vector<std::vector<double>> finite_diff_initial_position(
    const F& f, const double& t_in, const std::vector<double>& ts,
    const std::vector<double>& y_in, const std::vector<double>& theta,
    const std::vector<double>& x, const std::vector<int>& x_int,
    const size_t& param_index, const double& diff) {
  std::stringstream msgs;
  std::vector<double> y_in_ub(y_in.size());
  std::vector<double> y_in_lb(y_in.size());
  for (size_t i = 0; i < y_in.size(); i++) {
    if (i == param_index) {
      y_in_ub[i] = y_in[i] + diff;
      y_in_lb[i] = y_in[i] - diff;
    } else {
      y_in_ub[i] = y_in[i];
      y_in_lb[i] = y_in[i];
    }
  }

  std::vector<std::vector<double>> ode_res_ub;
  std::vector<std::vector<double>> ode_res_lb;

  ode_res_ub = stan::math::integrate_ode_rk45(f, y_in_ub, t_in, ts, theta, x,
                                              x_int, &msgs);
  ode_res_lb = stan::math::integrate_ode_rk45(f, y_in_lb, t_in, ts, theta, x,
                                              x_int, &msgs);

  std::vector<std::vector<double>> results(ts.size());

  for (size_t i = 0; i < ode_res_ub.size(); i++)
    for (size_t j = 0; j < ode_res_ub[j].size(); j++)
      results[i].push_back((ode_res_ub[i][j] - ode_res_lb[i][j]) / (2 * diff));
  return results;
}

// test integrate_ode_rk45 with initial positions as doubles and parameters
// as vars against finite differences
template <typename F>
void test_ode_finite_diff_dv(const F& f, const double& t_in,
                             const std::vector<double>& ts,
                             const std::vector<double>& y_in,
                             const std::vector<double>& theta,
                             const std::vector<double>& x,
                             const std::vector<int>& x_int, const double& diff,
                             const double& diff2) {
  std::stringstream msgs;

  std::vector<std::vector<std::vector<double>>> finite_diff_res(theta.size());
  for (size_t i = 0; i < theta.size(); i++)
    finite_diff_res[i]
        = finite_diff_params(f, t_in, ts, y_in, theta, x, x_int, i, diff);

  std::vector<double> grads_eff;

  std::vector<stan::math::var> theta_v;
  for (size_t i = 0; i < theta.size(); i++)
    theta_v.push_back(theta[i]);

  std::vector<std::vector<stan::math::var>> ode_res;

  ode_res = stan::math::integrate_ode_rk45(f, y_in, t_in, ts, theta_v, x, x_int,
                                           &msgs);

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = 0; j < y_in.size(); j++) {
      grads_eff.clear();
      ode_res[i][j].grad(theta_v, grads_eff);

      for (size_t k = 0; k < theta.size(); k++)
        EXPECT_NEAR(grads_eff[k], finite_diff_res[k][i][j], diff2)
            << "Gradient of integrate_ode_rk45 failed with initial positions"
            << " known and parameters unknown at time index " << i
            << ", equation index " << j << ", and parameter index: " << k;

      stan::math::set_zero_all_adjoints();
    }
  }
}

// test integrate_ode_rk45 with initial positions as vars and parameters
// as doubles against finite differences
template <typename F>
void test_ode_finite_diff_vd(const F& f, const double& t_in,
                             const std::vector<double>& ts,
                             const std::vector<double>& y_in,
                             const std::vector<double>& theta,
                             const std::vector<double>& x,
                             const std::vector<int>& x_int, const double& diff,
                             const double& diff2) {
  std::stringstream msgs;

  std::vector<std::vector<std::vector<double>>> finite_diff_res(y_in.size());
  for (size_t i = 0; i < y_in.size(); i++)
    finite_diff_res[i] = finite_diff_initial_position(f, t_in, ts, y_in, theta,
                                                      x, x_int, i, diff);

  std::vector<double> grads_eff;

  std::vector<stan::math::var> y_in_v;
  for (size_t k = 0; k < y_in.size(); k++)
    y_in_v.push_back(y_in[k]);

  std::vector<std::vector<stan::math::var>> ode_res;

  ode_res = stan::math::integrate_ode_rk45(f, y_in_v, t_in, ts, theta, x, x_int,
                                           &msgs);

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = 0; j < y_in.size(); j++) {
      grads_eff.clear();
      ode_res[i][j].grad(y_in_v, grads_eff);

      for (size_t k = 0; k < y_in.size(); k++)
        EXPECT_NEAR(grads_eff[k], finite_diff_res[k][i][j], diff2)
            << "Gradient of integrate_ode_rk45 failed with initial positions"
            << " unknown and parameters known at time index " << i
            << ", equation index " << j << ", and parameter index: " << k;

      stan::math::set_zero_all_adjoints();
    }
  }
}

// test integrate_ode_rk45 with initial positions as vars and parameters
// as vars against finite differences
template <typename F>
void test_ode_finite_diff_vv(const F& f, const double& t_in,
                             const std::vector<double>& ts,
                             const std::vector<double>& y_in,
                             const std::vector<double>& theta,
                             const std::vector<double>& x,
                             const std::vector<int>& x_int, const double& diff,
                             const double& diff2) {
  std::stringstream msgs;

  std::vector<std::vector<std::vector<double>>> finite_diff_res_y(y_in.size());
  for (size_t i = 0; i < y_in.size(); i++)
    finite_diff_res_y[i] = finite_diff_initial_position(
        f, t_in, ts, y_in, theta, x, x_int, i, diff);

  std::vector<std::vector<std::vector<double>>> finite_diff_res_p(theta.size());
  for (size_t i = 0; i < theta.size(); i++)
    finite_diff_res_p[i]
        = finite_diff_params(f, t_in, ts, y_in, theta, x, x_int, i, diff);

  std::vector<double> grads_eff;
  std::vector<stan::math::var> y_in_v;
  for (size_t i = 0; i < y_in.size(); i++)
    y_in_v.push_back(y_in[i]);

  std::vector<stan::math::var> vars = y_in_v;

  std::vector<stan::math::var> theta_v;
  for (size_t i = 0; i < theta.size(); i++)
    theta_v.push_back(theta[i]);

  for (size_t i = 0; i < theta_v.size(); i++)
    vars.push_back(theta_v[i]);

  std::vector<std::vector<stan::math::var>> ode_res;

  ode_res = stan::math::integrate_ode_rk45(f, y_in_v, t_in, ts, theta_v, x,
                                           x_int, &msgs);

  for (size_t i = 0; i < ts.size(); i++) {
    for (size_t j = 0; j < y_in.size(); j++) {
      grads_eff.clear();
      ode_res[i][j].grad(vars, grads_eff);

      for (size_t k = 0; k < theta.size(); k++)
        EXPECT_NEAR(grads_eff[k + y_in.size()], finite_diff_res_p[k][i][j],
                    diff2)
            << "Gradient of integrate_ode_rk45 failed with initial positions"
            << " unknown and parameters unknown for param at time index " << i
            << ", equation index " << j << ", and parameter index: " << k;
      for (size_t k = 0; k < y_in.size(); k++)
        EXPECT_NEAR(grads_eff[k], finite_diff_res_y[k][i][j], diff2)
            << "Gradient of integrate_ode_rk45 failed with initial positions"
            << " unknown and parameters known for initial position at time "
               "index "
            << i << ", equation index " << j << ", and parameter index: " << k;

      stan::math::set_zero_all_adjoints();
    }
  }
}

template <typename F, typename T1, typename T2>
void test_ode_error_conditions(F& f, const double& t0,
                               const std::vector<double>& ts,
                               const std::vector<T1>& y0,
                               const std::vector<T2>& theta,
                               const std::vector<double>& x,
                               const std::vector<int>& x_int) {
  using stan::math::integrate_ode_rk45;
  std::stringstream msgs;

  ASSERT_NO_THROW((integrate_ode_rk45(f, y0, t0, ts, theta, x, x_int, 0)));
  ASSERT_EQ("", msgs.str());

  msgs.clear();
  std::vector<T1> y0_bad;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0_bad, t0, ts, theta, x, x_int, &msgs)),
      std::invalid_argument, "initial state has size 0");
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  double t0_bad = ts[0] + 0.1;
  std::stringstream expected_msg;
  expected_msg << "initial time is " << t0_bad << ", but must be less than "
               << ts[0];
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0_bad, ts, theta, x, x_int, &msgs)),
      std::domain_error, expected_msg.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  std::vector<double> ts_bad;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts_bad, theta, x, x_int, &msgs)),
      std::invalid_argument, "times has size 0");
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  ts_bad.push_back(3);
  ts_bad.push_back(1);
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts_bad, theta, x, x_int, &msgs)),
      std::domain_error, "times is not a valid sorted vector");
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  std::vector<T2> theta_bad;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts, theta_bad, x, x_int, &msgs)),
      std::out_of_range, "vector");
  EXPECT_EQ("", msgs.str());

  if (x.size() > 0) {
    msgs.clear();
    std::vector<double> x_bad;
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(f, y0, t0, ts, theta, x_bad, x_int, &msgs)),
        std::out_of_range, "vector");
    EXPECT_EQ("", msgs.str());
  }

  if (x_int.size() > 0) {
    msgs.clear();
    std::vector<int> x_int_bad;
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(f, y0, t0, ts, theta, x, x_int_bad, &msgs)),
        std::out_of_range, "vector");
    EXPECT_EQ("", msgs.str());
  }
}

template <typename F, typename T1, typename T2>
void test_ode_error_conditions_nan(F& f, const double& t0,
                                   const std::vector<double>& ts,
                                   const std::vector<T1>& y0,
                                   const std::vector<T2>& theta,
                                   const std::vector<double>& x,
                                   const std::vector<int>& x_int) {
  using stan::math::integrate_ode_rk45;
  std::stringstream msgs;
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::stringstream expected_is_nan;
  expected_is_nan << "is " << nan;

  ASSERT_NO_THROW((integrate_ode_rk45(f, y0, t0, ts, theta, x, x_int, 0)));
  ASSERT_EQ("", msgs.str());

  msgs.clear();
  std::vector<T1> y0_bad = y0;
  y0_bad[0] = nan;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0_bad, t0, ts, theta, x, x_int, &msgs)),
      std::domain_error, "initial state");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0_bad, t0, ts, theta, x, x_int, &msgs)),
      std::domain_error, expected_is_nan.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  double t0_bad = nan;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0_bad, ts, theta, x, x_int, &msgs)),
      std::domain_error, "initial time");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0_bad, ts, theta, x, x_int, &msgs)),
      std::domain_error, expected_is_nan.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  std::vector<double> ts_bad = ts;
  ts_bad[0] = nan;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts_bad, theta, x, x_int, &msgs)),
      std::domain_error, "times");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts_bad, theta, x, x_int, &msgs)),
      std::domain_error, expected_is_nan.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  std::vector<T2> theta_bad = theta;
  theta_bad[0] = nan;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts, theta_bad, x, x_int, &msgs)),
      std::domain_error, "ode parameters and data");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts, theta_bad, x, x_int, &msgs)),
      std::domain_error, expected_is_nan.str());
  EXPECT_EQ("", msgs.str());

  if (x.size() > 0) {
    msgs.clear();
    std::vector<double> x_bad = x;
    x_bad[0] = nan;
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(f, y0, t0, ts, theta, x_bad, x_int, &msgs)),
        std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(f, y0, t0, ts, theta, x_bad, x_int, &msgs)),
        std::domain_error, expected_is_nan.str());
    EXPECT_EQ("", msgs.str());
  }
}

template <typename F, typename T1, typename T2>
void test_ode_error_conditions_inf(F& f, const double& t0,
                                   const std::vector<double>& ts,
                                   const std::vector<T1>& y0,
                                   const std::vector<T2>& theta,
                                   const std::vector<double>& x,
                                   const std::vector<int>& x_int) {
  using stan::math::integrate_ode_rk45;
  std::stringstream msgs;
  double inf = std::numeric_limits<double>::infinity();
  std::stringstream expected_is_inf;
  expected_is_inf << "is " << inf;
  std::stringstream expected_is_neg_inf;
  expected_is_neg_inf << "is " << -inf;

  ASSERT_NO_THROW((integrate_ode_rk45(f, y0, t0, ts, theta, x, x_int, 0)));
  ASSERT_EQ("", msgs.str());

  msgs.clear();
  std::vector<T1> y0_bad = y0;
  y0_bad[0] = inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0_bad, t0, ts, theta, x, x_int, &msgs)),
      std::domain_error, "initial state");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0_bad, t0, ts, theta, x, x_int, &msgs)),
      std::domain_error, expected_is_inf.str());
  y0_bad[0] = -inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0_bad, t0, ts, theta, x, x_int, &msgs)),
      std::domain_error, "initial state");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0_bad, t0, ts, theta, x, x_int, &msgs)),
      std::domain_error, expected_is_neg_inf.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  double t0_bad = inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0_bad, ts, theta, x, x_int, &msgs)),
      std::domain_error, "initial time");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0_bad, ts, theta, x, x_int, &msgs)),
      std::domain_error, expected_is_inf.str());
  t0_bad = -inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0_bad, ts, theta, x, x_int, &msgs)),
      std::domain_error, "initial time");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0_bad, ts, theta, x, x_int, &msgs)),
      std::domain_error, expected_is_neg_inf.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  std::vector<double> ts_bad = ts;
  ts_bad[0] = inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts_bad, theta, x, x_int, &msgs)),
      std::domain_error, "times");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts_bad, theta, x, x_int, &msgs)),
      std::domain_error, expected_is_inf.str());
  ts_bad[0] = -inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts_bad, theta, x, x_int, &msgs)),
      std::domain_error, "times");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts_bad, theta, x, x_int, &msgs)),
      std::domain_error, expected_is_neg_inf.str());
  EXPECT_EQ("", msgs.str());

  msgs.clear();
  std::vector<T2> theta_bad = theta;
  theta_bad[0] = inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts, theta_bad, x, x_int, &msgs)),
      std::domain_error, "ode parameters and data");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts, theta_bad, x, x_int, &msgs)),
      std::domain_error, expected_is_inf.str());
  theta_bad[0] = -inf;
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts, theta_bad, x, x_int, &msgs)),
      std::domain_error, "ode parameters and data");
  EXPECT_THROW_MSG(
      (integrate_ode_rk45(f, y0, t0, ts, theta_bad, x, x_int, &msgs)),
      std::domain_error, expected_is_neg_inf.str());
  EXPECT_EQ("", msgs.str());

  if (x.size() > 0) {
    msgs.clear();
    std::vector<double> x_bad = x;
    x_bad[0] = inf;
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(f, y0, t0, ts, theta, x_bad, x_int, &msgs)),
        std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(f, y0, t0, ts, theta, x_bad, x_int, &msgs)),
        std::domain_error, expected_is_inf.str());
    x_bad[0] = -inf;
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(f, y0, t0, ts, theta, x_bad, x_int, &msgs)),
        std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG(
        (integrate_ode_rk45(f, y0, t0, ts, theta, x_bad, x_int, &msgs)),
        std::domain_error, expected_is_neg_inf.str());
    EXPECT_EQ("", msgs.str());
  }
}

template <typename F>
void test_ode_error_conditions_vd(const F& f, const double& t_in,
                                  const std::vector<double>& ts,
                                  const std::vector<double>& y_in,
                                  const std::vector<double>& theta,
                                  const std::vector<double>& x,
                                  const std::vector<int>& x_int) {
  std::vector<stan::math::var> y_var;
  for (size_t i = 0; i < y_in.size(); i++)
    y_var.push_back(y_in[i]);
  test_ode_error_conditions(f, t_in, ts, y_var, theta, x, x_int);
  test_ode_error_conditions_nan(f, t_in, ts, y_var, theta, x, x_int);
}
template <typename F>
void test_ode_error_conditions_dv(const F& f, const double& t_in,
                                  const std::vector<double>& ts,
                                  const std::vector<double>& y_in,
                                  const std::vector<double>& theta,
                                  const std::vector<double>& x,
                                  const std::vector<int>& x_int) {
  std::vector<stan::math::var> theta_var;
  for (size_t i = 0; i < theta.size(); i++)
    theta_var.push_back(theta[i]);
  test_ode_error_conditions(f, t_in, ts, y_in, theta_var, x, x_int);
  test_ode_error_conditions_nan(f, t_in, ts, y_in, theta_var, x, x_int);
}
template <typename F>
void test_ode_error_conditions_vv(const F& f, const double& t_in,
                                  const std::vector<double>& ts,
                                  const std::vector<double>& y_in,
                                  const std::vector<double>& theta,
                                  const std::vector<double>& x,
                                  const std::vector<int>& x_int) {
  std::vector<stan::math::var> y_var;
  for (size_t i = 0; i < y_in.size(); i++)
    y_var.push_back(y_in[i]);
  std::vector<stan::math::var> theta_var;
  for (size_t i = 0; i < theta.size(); i++)
    theta_var.push_back(theta[i]);
  test_ode_error_conditions(f, t_in, ts, y_var, theta_var, x, x_int);
  test_ode_error_conditions_nan(f, t_in, ts, y_var, theta_var, x, x_int);
  test_ode_error_conditions_inf(f, t_in, ts, y_var, theta_var, x, x_int);
}

template <typename F>
void test_ode(const F& f, const double& t_in, const std::vector<double>& ts,
              const std::vector<double>& y_in, const std::vector<double>& theta,
              const std::vector<double>& x, const std::vector<int>& x_int,
              const double& diff, const double& diff2) {
  test_ode_finite_diff_vd(f, t_in, ts, y_in, theta, x, x_int, diff, diff2);
  test_ode_finite_diff_dv(f, t_in, ts, y_in, theta, x, x_int, diff, diff2);
  test_ode_finite_diff_vv(f, t_in, ts, y_in, theta, x, x_int, diff, diff2);

  test_ode_error_conditions_vd(f, t_in, ts, y_in, theta, x, x_int);
  test_ode_error_conditions_dv(f, t_in, ts, y_in, theta, x, x_int);
  test_ode_error_conditions_vv(f, t_in, ts, y_in, theta, x, x_int);
}
