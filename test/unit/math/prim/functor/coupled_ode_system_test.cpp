#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/functor/mock_ode_functor.hpp>
#include <test/unit/math/prim/functor/mock_throwing_ode_functor.hpp>
#include <vector>
#include <string>

struct StanMathCoupledOdeSystem : public ::testing::Test {
  std::stringstream msgs;
  std::vector<double> x;
  std::vector<int> x_int;
};

TEST_F(StanMathCoupledOdeSystem, initial_state_dd) {
  using stan::math::coupled_ode_system;
  mock_ode_functor base_ode;

  const int N = 3;
  const int M = 4;

  Eigen::VectorXd y0_d(N);
  std::vector<double> theta_d(M);

  for (int n = 0; n < N; n++)
    y0_d(n) = n + 1;
  for (int m = 0; m < M; m++)
    theta_d[m] = 10 * (m + 1);

  coupled_ode_system<mock_ode_functor, double, std::vector<double>,
                     std::vector<double>, std::vector<int>>
      coupled_system_dd(base_ode, y0_d, &msgs, theta_d, x, x_int);

  std::vector<double> state = coupled_system_dd.initial_state();
  for (int n = 0; n < N; n++)
    EXPECT_FLOAT_EQ(y0_d(n), state[n])
        << "we don't need derivatives of y0; "
        << "initial state gets the initial values";
  for (size_t n = N; n < state.size(); n++)
    EXPECT_FLOAT_EQ(0.0, state[n]);
}

TEST_F(StanMathCoupledOdeSystem, size) {
  using stan::math::coupled_ode_system;
  mock_ode_functor base_ode;

  const int N = 3;
  const int M = 4;

  Eigen::VectorXd y0_d(N);
  std::vector<double> theta_d(M, 0.0);

  coupled_ode_system<mock_ode_functor, double, int, double, Eigen::MatrixXd>
      coupled_system_dd(base_ode, y0_d, &msgs, 1, 1.0, y0_d);

  EXPECT_EQ(N, coupled_system_dd.size());
}

TEST_F(StanMathCoupledOdeSystem, recover_exception) {
  using stan::math::coupled_ode_system;
  std::string message = "ode throws";

  const int N = 3;
  const int M = 4;

  mock_throwing_ode_functor<std::logic_error> throwing_ode(message);

  Eigen::VectorXd y0_d(N);
  std::vector<double> theta_v(M);

  coupled_ode_system<mock_throwing_ode_functor<std::logic_error>, double,
                     std::vector<double>, std::vector<double>, std::vector<int>>
      coupled_system_dd(throwing_ode, y0_d, &msgs, theta_v, x, x_int);

  std::vector<double> y(3);
  std::vector<double> dy_dt(3);

  double t = 10;

  EXPECT_THROW_MSG(coupled_system_dd(y, dy_dt, t), std::logic_error, message);
}
