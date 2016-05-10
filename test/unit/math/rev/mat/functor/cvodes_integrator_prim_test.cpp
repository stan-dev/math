#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/arr/functor/mock_ode_functor.hpp>
#include <test/unit/math/prim/arr/functor/mock_throwing_ode_functor.hpp>

struct StanMathOdeCVode : public ::testing::Test {
  void SetUp() {
    stan::math::recover_memory();
    ts = std::vector<double>(1,10);
    t0 = 0;
  }
  std::stringstream msgs;
  std::vector<double> x;
  std::vector<int> x_int;
  double t0;
  std::vector<double> ts;
};
TEST_F(StanMathOdeCVode, recover_exception) {
  using stan::math::cvodes_integrator;
  std::string message = "ode throws";

  const int N = 3;
  const int M = 4;

  mock_throwing_ode_functor<std::logic_error> throwing_ode(message);

  std::vector<double> y0_d(N, 0.0);
  std::vector<double> theta_v(M, 0.0);

  cvodes_integrator<mock_throwing_ode_functor<std::logic_error>, double, double>
    integrator_dd(throwing_ode, y0_d, t0, theta_v, x, x_int, ts, 1e-8, 1e-10, 1e6, 1, &msgs);

  std::vector<double> y(3,0);
  std::vector<double> dy_dt(3,0);

  N_Vector nv_y = N_VMake_Serial(N, &y[0]);
  N_Vector nv_dy_dt = N_VMake_Serial(N, &dy_dt[0]);


  double t = 10;

  EXPECT_THROW_MSG(integrator_dd.ode_rhs(t, nv_y, nv_dy_dt, &integrator_dd),
                   std::logic_error,
                   message);
}
