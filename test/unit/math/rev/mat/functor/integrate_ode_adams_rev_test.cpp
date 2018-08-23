#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
#include <test/unit/math/rev/mat/functor/util_cvodes_adams.hpp>
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

  std::vector<std::vector<var> > ode_res_vd = stan::math::integrate_ode_adams(
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

  EXPECT_THROW_MSG(stan::math::integrate_ode_adams(
                       harm_osc, promote_scalar<T_y0>(y0), t0, ts,
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

// TODO(Yi Zhang): failure
// TEST(StanAgradRevOde_integrate_ode, lorenz_finite_diff) {
//   lorenz_ode_fun lorenz;

//   std::vector<double> y0;
//   std::vector<double> theta;
//   double t0;
//   std::vector<double> ts;

//   t0 = 0;

//   theta.push_back(10.0);
//   theta.push_back(28.0);
//   theta.push_back(8.0 / 3.0);
//   y0.push_back(10.0);
//   y0.push_back(1.0);
//   y0.push_back(1.0);

//   std::vector<double> x;
//   std::vector<int> x_int;

//   for (int i = 0; i < 100; i++)
//     ts.push_back(0.1 * (i + 1));

//   test_ode_cvode(lorenz, t0, ts, y0, theta, x, x_int, 1e-8, 1e-1);
// }
