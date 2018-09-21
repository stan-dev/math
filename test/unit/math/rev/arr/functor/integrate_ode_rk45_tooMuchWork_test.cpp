#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>
#include <boost/numeric/odeint.hpp>
// very small michaelis menten example
#include <test/unit/math/rev/arr/functor/coupled_mm.hpp>
#include <test/unit/util.hpp>
#include <vector>

// test which triggers the too much work exception from odeint
TEST(StanOde_tooMuchWork_test, odeint_coupled_mm) {
  coupled_mm_ode_fun f_;

  // initial value and parameters from model definition
  std::vector<stan::math::var> y0_v(2);
  y0_v[0] = 1E5;
  y0_v[1] = 1E-1;

  double t0 = 0;

  std::vector<double> ts_long;
  ts_long.push_back(1E10);

  std::vector<double> ts_short;
  ts_short.push_back(1);

  std::vector<stan::math::var> theta_v(4);

  theta_v[0] = 1.0;
  theta_v[1] = 0.5;
  theta_v[2] = 0.5;
  theta_v[3] = 0.1;

  std::vector<double> data;

  std::vector<int> data_int;

  EXPECT_THROW_MSG(
      stan::math::integrate_ode_rk45(f_, y0_v, t0, ts_long, theta_v, data,
                                     data_int, 0, 1E-6, 1E-6, 100),
      boost::numeric::odeint::no_progress_error,
      "Max number of iterations exceeded (100).");

  EXPECT_NO_THROW(stan::math::integrate_ode_rk45(
      f_, y0_v, t0, ts_short, theta_v, data, data_int, 0, 1E-6, 1E-6, 100));
}
