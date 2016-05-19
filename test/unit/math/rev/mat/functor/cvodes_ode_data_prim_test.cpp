#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/arr/functor/mock_ode_functor.hpp>
#include <test/unit/math/prim/arr/functor/mock_throwing_ode_functor.hpp>

struct StanMathOdeCVode : public ::testing::Test {
  void SetUp() {
    stan::math::recover_memory();
    t0 = 0;
  }
  std::stringstream msgs;
  std::vector<double> x;
  std::vector<int> x_int;
  double t0;
};
TEST_F(StanMathOdeCVode, recover_exception) {
  using stan::math::cvodes_ode_data;
  std::string message = "ode throws";

  const int N = 3;
  const int M = 4;

  mock_throwing_ode_functor<std::logic_error> throwing_ode(message);

  std::vector<double> y0_d(N, 0.0);
  std::vector<double> theta_v(M, 0.0);

  cvodes_ode_data<mock_throwing_ode_functor<std::logic_error>, double, double>
    ode_data_dd(throwing_ode, y0_d, theta_v, x, x_int, &msgs);

  std::vector<double> y(3,0);
  std::vector<double> dy_dt(3,0);

  N_Vector nv_y = N_VMake_Serial(N, &y[0]);
  N_Vector nv_dy_dt = N_VMake_Serial(N, &dy_dt[0]);


  double t = 10;

  EXPECT_THROW_MSG(ode_data_dd.ode_rhs(t, nv_y, nv_dy_dt, &ode_data_dd),
                   std::logic_error,
                   message);
}
