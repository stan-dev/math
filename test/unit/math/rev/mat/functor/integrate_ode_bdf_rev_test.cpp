#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <test/unit/math/rev/mat/functor/util_cvodes_bdf.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/arr/functor/lorenz.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

template <typename F, typename T_y0, typename T_theta>
void sho_value_test(F harm_osc, std::vector<double>& y0, double t0,
                    std::vector<double>& ts, std::vector<double>& theta,
                    std::vector<double>& x, std::vector<int>& x_int) {
  using stan::math::promote_scalar;
  using stan::math::var;

  std::vector<std::vector<var> > ode_res_vd = stan::math::integrate_ode_bdf(
      harm_osc, promote_scalar<T_y0>(y0), t0, ts,
      promote_scalar<T_theta>(theta), x, x_int);

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

  test_ode_cvode(harm_osc, t0, ts, y0, theta, x, x_int, 1e-8, 1e-4);

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

  test_ode_cvode(harm_osc, t0, ts, y0, theta, x, x_int, 1e-8, 1e-4);

  sho_value_test<harm_osc_ode_data_fun, double, var>(harm_osc, y0, t0, ts,
                                                     theta, x, x_int);
  sho_value_test<harm_osc_ode_data_fun, var, double>(harm_osc, y0, t0, ts,
                                                     theta, x, x_int);
  sho_value_test<harm_osc_ode_data_fun, var, var>(harm_osc, y0, t0, ts, theta,
                                                  x, x_int);
}

template <typename T_y0, typename T_theta, typename F>
void sho_error_test(F harm_osc, std::vector<double>& y0, double t0,
                    std::vector<double>& ts, std::vector<double>& theta,
                    std::vector<double>& x, std::vector<int>& x_int,
                    std::string error_msg) {
  using stan::math::promote_scalar;
  using stan::math::var;

  EXPECT_THROW_MSG(
      stan::math::integrate_ode_bdf(harm_osc, promote_scalar<T_y0>(y0), t0, ts,
                                    promote_scalar<T_theta>(theta), x, x_int),
      std::invalid_argument, error_msg);
}

// TODO(carpenter): g++6 failure
TEST(StanAgradRevOde_integrate_ode, harmonic_oscillator_finite_diff) {
  sho_finite_diff_test(0);
  sho_finite_diff_test(2.0);
  sho_finite_diff_test(-2.0);

  sho_data_finite_diff_test(0);
  sho_data_finite_diff_test(2.5);
  sho_data_finite_diff_test(-2.5);
}

TEST(StanAgradRevOde_integrate_ode, harmonic_oscillator_error) {
  using stan::math::var;
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

  // aligned error handling with non-stiff case
  std::string error_msg
      = "cvodes_ode_data: dz_dt (3) and states (2) must match in size";

  sho_error_test<double, var>(harm_osc, y0, t0, ts, theta, x, x_int, error_msg);
  sho_error_test<var, double>(harm_osc, y0, t0, ts, theta, x, x_int, error_msg);
  sho_error_test<var, var>(harm_osc, y0, t0, ts, theta, x, x_int, error_msg);
}

// TODO(carpenter): g++6 failure
TEST(StanAgradRevOde_integrate_ode, lorenz_finite_diff) {
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

  test_ode_cvode(lorenz, t0, ts, y0, theta, x, x_int, 1e-8, 1e-1);
}

TEST(StanAgradRevOde_integrate_ode_bdf, time_steps_as_param) {
  using stan::math::integrate_ode_bdf;
  using stan::math::to_var;

  const double t0 = 0.0;
  harm_osc_ode_fun ode;
  std::vector<double> theta{0.15};
  std::vector<double> y0{1.0, 0.0};
  std::vector<stan::math::var> ts;
  for (int i = 0; i < 100; i++)
    ts.push_back(t0 + 0.1 * (i + 1));
  std::vector<double> x;
  std::vector<int> x_int;
  std::vector<stan::math::var> y0v = to_var(y0);
  std::vector<stan::math::var> thetav = to_var(theta);

  std::vector<std::vector<stan::math::var> > res;
  auto test_val = [&res]() {
    EXPECT_NEAR(0.995029, res[0][0].val(), 1e-5);
    EXPECT_NEAR(-0.0990884, res[0][1].val(), 1e-5);
    EXPECT_NEAR(-0.421907, res[99][0].val(), 1e-5);
    EXPECT_NEAR(0.246407, res[99][1].val(), 1e-5);
  };
  res = integrate_ode_bdf(ode, y0, t0, ts, theta, x, x_int);
  test_val();
  res = integrate_ode_bdf(ode, y0v, t0, ts, theta, x, x_int);
  test_val();
  res = integrate_ode_bdf(ode, y0, t0, ts, thetav, x, x_int);
  test_val();
  res = integrate_ode_bdf(ode, y0v, t0, ts, thetav, x, x_int);
  test_val();
}

TEST(StanAgradRevOde_integrate_ode_bdf, time_steps_as_param_AD) {
  using stan::math::integrate_ode_bdf;
  using stan::math::to_var;
  using stan::math::value_of;
  using stan::math::var;
  const double t0 = 0.0;
  const int nt = 100;  // nb. of time steps
  const int ns = 2;    // nb. of states
  std::ostream* msgs = NULL;

  harm_osc_ode_fun ode;

  std::vector<double> theta{0.15};
  std::vector<double> y0{1.0, 0.0};
  std::vector<stan::math::var> ts;
  for (int i = 0; i < nt; i++)
    ts.push_back(t0 + 0.1 * (i + 1));

  std::vector<double> x;
  std::vector<int> x_int;
  std::vector<stan::math::var> y0v = to_var(y0);
  std::vector<stan::math::var> thetav = to_var(theta);

  std::vector<std::vector<stan::math::var> > res;
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
  res = integrate_ode_bdf(ode, y0, t0, ts, theta, x, x_int);
  test_ad();
  res = integrate_ode_bdf(ode, y0v, t0, ts, theta, x, x_int);
  test_ad();
  res = integrate_ode_bdf(ode, y0, t0, ts, thetav, x, x_int);
  test_ad();
  res = integrate_ode_bdf(ode, y0v, t0, ts, thetav, x, x_int);
  test_ad();
}
