#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <test/unit/math/rev/functor/util_rk45.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/functor/forced_harmonic_oscillator.hpp>
#include <test/unit/math/prim/functor/lorenz.hpp>
#include <iostream>
#include <sstream>
#include <vector>

template <typename F, typename T_y0, typename T_theta>
void sho_value_test(F harm_osc, std::vector<double>& y0, double t0,
                    std::vector<double>& ts, std::vector<double>& theta,
                    std::vector<double>& x, std::vector<int>& x_int) {
  using stan::math::promote_scalar;
  using stan::math::var;

  std::vector<std::vector<var>> ode_res_vd = stan::math::integrate_ode_rk45(
      harm_osc, promote_scalar<T_y0>(y0), t0, ts,
      promote_scalar<T_theta>(theta), x, x_int, 0);
  EXPECT_NEAR(0.995029, ode_res_vd[0][0].val(), 1e-5);
  EXPECT_NEAR(-0.0990884, ode_res_vd[0][1].val(), 1e-5);

  EXPECT_NEAR(-0.421907, ode_res_vd[99][0].val(), 1e-5);
  EXPECT_NEAR(0.246407, ode_res_vd[99][1].val(), 1e-5);
}

void sho_finite_diff_test(double t0) {
  using stan::math::var;
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

  test_ode(harm_osc, t0, ts, y0, theta, x, x_int, 1e-8, 1e-4);

  sho_value_test<harm_osc_ode_fun, double, var>(harm_osc, y0, t0, ts, theta, x,
                                                x_int);
  sho_value_test<harm_osc_ode_fun, var, double>(harm_osc, y0, t0, ts, theta, x,
                                                x_int);
  sho_value_test<harm_osc_ode_fun, var, var>(harm_osc, y0, t0, ts, theta, x,
                                             x_int);
}

void sho_data_finite_diff_test(double t0) {
  using stan::math::var;
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

  test_ode(harm_osc, t0, ts, y0, theta, x, x_int, 1e-8, 1e-4);

  sho_value_test<harm_osc_ode_data_fun, double, var>(harm_osc, y0, t0, ts,
                                                     theta, x, x_int);
  sho_value_test<harm_osc_ode_data_fun, var, double>(harm_osc, y0, t0, ts,
                                                     theta, x, x_int);
  sho_value_test<harm_osc_ode_data_fun, var, var>(harm_osc, y0, t0, ts, theta,
                                                  x, x_int);
}

TEST(StanAgradRevOde_integrate_ode_rk45, harmonic_oscillator_finite_diff) {
  sho_finite_diff_test(0);
  sho_finite_diff_test(1.0);
  sho_finite_diff_test(-1.0);

  sho_data_finite_diff_test(0);
  sho_data_finite_diff_test(1.0);
  sho_data_finite_diff_test(-1.0);
}

TEST(StanAgradRevOde_integrate_ode_rk45, lorenz_finite_diff) {
  lorenz_ode_fun lorenz;

  std::vector<double> y0;
  std::vector<double> theta;
  double t0;
  std::vector<double> ts;

  t0 = 0;

  theta.push_back(10.0);
  theta.push_back(28.0);
  theta.push_back(8.0 / 3.0);
  y0.push_back(10.0);
  y0.push_back(1.0);
  y0.push_back(1.0);

  std::vector<double> x;
  std::vector<int> x_int;

  for (int i = 0; i < 100; i++)
    ts.push_back(0.1 * (i + 1));

  test_ode(lorenz, t0, ts, y0, theta, x, x_int, 1e-8, 1e-1);
}

TEST(StanAgradRevOde_integrate_ode_rk45, time_steps_as_param) {
  using stan::math::integrate_ode_rk45;
  using stan::math::to_var;
  using stan::math::value_of;

  const double t0 = 0.0;
  forced_harm_osc_ode_fun ode;
  std::vector<double> theta{0.15, 0.25};
  std::vector<double> y0{1.0, 0.0};
  std::vector<stan::math::var> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));
  std::vector<double> x;
  std::vector<int> x_int;
  std::vector<stan::math::var> y0v = to_var(y0);
  std::vector<stan::math::var> thetav = to_var(theta);
  stan::math::var t0v = 0.0;

  std::vector<std::vector<stan::math::var>> res;

  std::vector<std::vector<double>> res_d
      = integrate_ode_rk45(ode, y0, t0, value_of(ts), theta, x, x_int);

  // here we only test first & last steps, and rely on the
  // fact that results in-between affect the initial
  // condition of the last step to check their validity.
  auto test_val = [&res_d, &res]() {
    EXPECT_NEAR(res_d[0][0], res[0][0].val(), 1e-5);
    EXPECT_NEAR(res_d[0][1], res[0][1].val(), 1e-5);
    EXPECT_NEAR(res_d[99][0], res[99][0].val(), 1e-5);
    EXPECT_NEAR(res_d[99][1], res[99][1].val(), 1e-5);
  };
  res = integrate_ode_rk45(ode, y0, t0, ts, theta, x, x_int);
  test_val();
  res = integrate_ode_rk45(ode, y0v, t0, ts, theta, x, x_int);
  test_val();
  res = integrate_ode_rk45(ode, y0, t0, ts, thetav, x, x_int);
  test_val();
  res = integrate_ode_rk45(ode, y0v, t0, ts, thetav, x, x_int);
  test_val();
  res = integrate_ode_rk45(ode, y0, t0v, ts, theta, x, x_int);
  test_val();
  res = integrate_ode_rk45(ode, y0v, t0v, ts, theta, x, x_int);
  test_val();
  res = integrate_ode_rk45(ode, y0, t0v, ts, thetav, x, x_int);
  test_val();
  res = integrate_ode_rk45(ode, y0v, t0v, ts, thetav, x, x_int);
  test_val();
}

TEST(StanAgradRevOde_integrate_ode_rk45, time_steps_as_param_AD) {
  using stan::math::integrate_ode_rk45;
  using stan::math::to_var;
  using stan::math::value_of;
  using stan::math::var;
  const double t0 = 0.0;
  const int nt = 100;  // nb. of time steps
  const int ns = 2;    // nb. of states
  std::ostream* msgs = NULL;

  forced_harm_osc_ode_fun ode;

  std::vector<double> theta{0.15, 0.25};
  std::vector<double> y0{1.0, 0.0};
  std::vector<stan::math::var> ts;
  for (int i = 0; i < nt; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x;
  std::vector<int> x_int;
  std::vector<stan::math::var> y0v = to_var(y0);
  std::vector<stan::math::var> thetav = to_var(theta);
  stan::math::var t0v = 0.0;

  std::vector<std::vector<stan::math::var>> res;
  std::vector<double> g;
  auto test_ad = [&res, &g, &ts, &ode, &nt, &ns, &theta, &x, &x_int, &msgs]() {
    for (auto i = 0; i < nt; ++i) {
      std::vector<double> res_d = value_of(res[i]);
      for (auto j = 0; j < ns; ++j) {
        g.clear();
        res[i][j].grad(ts, g);
        for (auto k = 0; k < nt; ++k) {
          if (k != i) {
            EXPECT_FLOAT_EQ(g[k], 0.0);
          } else {
            std::vector<double> y0(res_d.begin(), res_d.begin() + ns);
            EXPECT_FLOAT_EQ(g[k],
                            ode(ts[i].val(), y0, theta, x, x_int, msgs)[j]);
          }
        }
        stan::math::set_zero_all_adjoints();
      }
    }
  };
  res = integrate_ode_rk45(ode, y0, t0, ts, theta, x, x_int);
  test_ad();
  res = integrate_ode_rk45(ode, y0v, t0, ts, theta, x, x_int);
  test_ad();
  res = integrate_ode_rk45(ode, y0, t0, ts, thetav, x, x_int);
  test_ad();
  res = integrate_ode_rk45(ode, y0v, t0, ts, thetav, x, x_int);
  test_ad();
}

TEST(StanAgradRevOde_integrate_ode_rk45, t0_as_param_AD) {
  using stan::math::integrate_ode_rk45;
  using stan::math::to_var;
  using stan::math::value_of;
  using stan::math::var;
  const double t0 = 0.0;
  std::ostream* msgs = NULL;

  harm_osc_ode_fun ode;

  std::vector<double> theta{0.15};
  std::vector<double> y0{1.0, 0.0};
  std::vector<double> ts = {5.0, 10.0};

  std::vector<double> x;
  std::vector<int> x_int;
  std::vector<stan::math::var> y0v = to_var(y0);
  std::vector<stan::math::var> thetav = to_var(theta);
  stan::math::var t0v = to_var(t0);

  std::vector<std::vector<stan::math::var>> res;
  auto test_ad = [&res, &t0v, &ode, &theta, &x, &x_int, &msgs]() {
    res[0][0].grad();
    EXPECT_FLOAT_EQ(t0v.adj(), -0.66360742442816977871);
    stan::math::set_zero_all_adjoints();
    res[0][1].grad();
    EXPECT_FLOAT_EQ(t0v.adj(), 0.23542843380353062344);
    stan::math::set_zero_all_adjoints();
    res[1][0].grad();
    EXPECT_FLOAT_EQ(t0v.adj(), -0.2464078910913158893);
    stan::math::set_zero_all_adjoints();
    res[1][1].grad();
    EXPECT_FLOAT_EQ(t0v.adj(), -0.38494826636037426937);
    stan::math::set_zero_all_adjoints();
  };
  res = integrate_ode_rk45(ode, y0, t0v, ts, theta, x, x_int, nullptr, 1e-10,
                           1e-10, 1e6);
  test_ad();
  res = integrate_ode_rk45(ode, y0v, t0v, ts, theta, x, x_int, nullptr, 1e-10,
                           1e-10, 1e6);
  test_ad();
  res = integrate_ode_rk45(ode, y0, t0v, ts, thetav, x, x_int, nullptr, 1e-10,
                           1e-10, 1e6);
  test_ad();
  res = integrate_ode_rk45(ode, y0v, t0v, ts, thetav, x, x_int, nullptr, 1e-10,
                           1e-10, 1e6);
  test_ad();
}
