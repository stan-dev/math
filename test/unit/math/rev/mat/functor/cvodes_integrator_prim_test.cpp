#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/arr/functor/mock_ode_functor.hpp>
#include <test/unit/math/prim/arr/functor/mock_throwing_ode_functor.hpp>

struct StanMathOdeCVode : public ::testing::Test {
  std::stringstream msgs;
  std::vector<double> x;
  std::vector<int> x_int;
  double t0;
};

TEST_F(StanMathOdeCVode, decouple_ode_states_dd) {
  using stan::math::cvodes_integrator;
  using stan::math::decouple_ode_states;

  std::vector<double> y0(2);
  y0[0] = 1.0;
  y0[1] = 0.5;

  std::vector<double> theta(1);
  theta[0] = 0.15;

  stan::math::ode_model<harm_osc_ode_fun> harm_ode(harm_osc_ode_fun(), theta, x, x_int, &msgs);

  cvodes_integrator<harm_osc_ode_fun>
    integrator(harm_ode, y0, t0, false, false, 1e-8, 1e-10, 1e6, 1);

  int T = 10;
  int k = 0;
  std::vector<std::vector<double> > ys_coupled(T);
  for (int t = 0; t < T; t++) {
    std::vector<double> coupled_state(integrator.size(), 0.0);
    for (int n = 0; n < integrator.size(); n++)
      coupled_state[n] = ++k;
    ys_coupled[t] = coupled_state;
  }

  std::vector<std::vector<double> > ys;
  ys = decouple_ode_states(ys_coupled, y0, theta);

  ASSERT_EQ(T, ys.size());
  for (int t = 0; t < T; t++)
    ASSERT_EQ(2, ys[t].size());

  for (int t = 0; t < T; t++)
    for (int n = 0; n < 2; n++)
      EXPECT_FLOAT_EQ(ys_coupled[t][n], ys[t][n])
        << "(" << n << "," << t << "): "
        << "for (double, double) the coupled system is the base system";
}
/* obsolete
TEST_F(StanMathOdeCVode, initial_state_dd) {
  using stan::math::cvodes_integrator;

  const int N = 3;
  const int M = 4;

  std::vector<double> y0_d(N, 0.0);
  std::vector<double> theta_d(M, 0.0);

  for (int n = 0; n < N; n++)
    y0_d[n] = n+1;
  for (int m = 0; m < M; m++)
    theta_d[m] = 10 * (m+1);

  stan::math::ode_model<mock_ode_functor> mock_ode(mock_ode_functor(), theta_d, x, x_int, &msgs);
  
  cvodes_integrator<mock_ode_functor>
    integrator_dd(mock_ode, y0_d, t0, false, false, 1e-8, 1e-10, 1e6, 1);

  std::vector<double> state  = integrator_dd.initial_state();
  for (int n = 0; n < N; n++)
    EXPECT_FLOAT_EQ(y0_d[n], state[n])
      << "initial state gets the initial values";
  for (size_t n = N; n < state.size(); n++)
    EXPECT_FLOAT_EQ(0.0, state[n]);
}
*/
TEST_F(StanMathOdeCVode, size) {
  using stan::math::cvodes_integrator;

  const int N = 3;
  const int M = 4;

  std::vector<double> y0_d(N, 0.0);
  std::vector<double> theta_d(M, 0.0);

  stan::math::ode_model<mock_ode_functor> mock_ode(mock_ode_functor(), theta_d, x, x_int, &msgs);
  
  cvodes_integrator<mock_ode_functor>
    integrator_dd(mock_ode, y0_d, t0, false, false, 1e-8, 1e-10, 1e6, 1);

  EXPECT_EQ(N, integrator_dd.size());
}


TEST_F(StanMathOdeCVode, recover_exception) {
  using stan::math::cvodes_integrator;
  std::string message = "ode throws";

  const int N = 3;
  const int M = 4;

  mock_throwing_ode_functor<std::logic_error> throwing_ode(message);

  std::vector<double> y0_d(N, 0.0);
  std::vector<double> theta_d(M, 0.0);

  typedef mock_throwing_ode_functor<std::logic_error> throw_mock_t;
  stan::math::ode_model<throw_mock_t> mock_ode(throw_mock_t(message, 1), theta_d, x, x_int, &msgs);
  
  cvodes_integrator<throw_mock_t>
    integrator_dd(mock_ode, y0_d, t0, false, false, 1e-8, 1e-10, 1e6, 1);

  std::vector<double> y(3,0);
  std::vector<double> dy_dt(3,0);

  double t = 10;

  EXPECT_THROW_MSG(integrator_dd.rhs(&y[0], &dy_dt[0], t),
                   std::logic_error,
                   message);
}
