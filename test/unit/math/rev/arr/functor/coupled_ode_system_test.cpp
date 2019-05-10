#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/arr/functor/mock_ode_functor.hpp>
#include <test/unit/math/prim/arr/functor/mock_throwing_ode_functor.hpp>
#include <vector>
#include <string>

struct StanAgradRevOde : public ::testing::Test {
  void SetUp() { stan::math::recover_memory(); }
  std::stringstream msgs;
  std::vector<double> x;
  std::vector<int> x_int;
};

// ******************** DV ****************************
TEST_F(StanAgradRevOde, coupled_ode_system_dv) {
  using stan::math::coupled_ode_system;

  stan::math::start_nested();

  harm_osc_ode_fun harm_osc;

  std::vector<stan::math::var> theta;
  std::vector<double> z0;
  std::vector<double> y0;
  double t0;
  const size_t z_size = 4;
  std::vector<double> dz_dt(z_size, 0);

  double gamma(0.15);
  t0 = 0;

  theta.push_back(gamma);
  y0.push_back(1.0);
  y0.push_back(0.5);

  z0.push_back(1.0);
  z0.push_back(0.5);
  z0.push_back(1.0);
  z0.push_back(2.0);

  std::size_t stack_size = stan::math::nested_size();

  coupled_ode_system<harm_osc_ode_fun, double, stan::math::var> system(
      harm_osc, y0, theta, x, x_int, &msgs);

  EXPECT_EQ(stack_size, stan::math::nested_size())
      << "expecting no new things on the stack";

  system(z0, dz_dt, t0);

  EXPECT_FLOAT_EQ(0.5, dz_dt[0]);
  EXPECT_FLOAT_EQ(-1.075, dz_dt[1]);
  EXPECT_FLOAT_EQ(2, dz_dt[2]);
  EXPECT_FLOAT_EQ(-1.8, dz_dt[3]);

  stan::math::recover_memory_nested();
}
TEST_F(StanAgradRevOde, initial_state_dv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * M;

  std::vector<double> y0_d(N, 0.0);
  std::vector<var> theta_v(M, 0.0);

  for (size_t n = 0; n < N; n++)
    y0_d[n] = n + 1;
  for (size_t m = 0; m < M; m++)
    theta_v[m] = 10 * (m + 1);

  coupled_ode_system<mock_ode_functor, double, var> coupled_system_dv(
      base_ode, y0_d, theta_v, x, x_int, &msgs);

  std::vector<double> state = coupled_system_dv.initial_state();
  for (size_t n = 0; n < N; n++)
    EXPECT_FLOAT_EQ(y0_d[n], state[n])
        << "we don't need derivatives of y0; "
        << "initial state gets the initial values";
  for (size_t n = N; n < state.size(); n++)
    EXPECT_FLOAT_EQ(0.0, state[n]);
}
TEST_F(StanAgradRevOde, size_dv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * M;

  std::vector<double> y0_d(N, 0.0);
  std::vector<var> theta_v(M, 0.0);

  coupled_ode_system<mock_ode_functor, double, var> coupled_system_dv(
      base_ode, y0_d, theta_v, x, x_int, &msgs);

  EXPECT_EQ(z_size, coupled_system_dv.size());
}

TEST_F(StanAgradRevOde, memory_recovery_dv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * M;

  std::vector<double> y0_d(N, 0.0);
  std::vector<var> theta_v(M, 0.0);

  coupled_ode_system<mock_ode_functor, double, var> coupled_system_dv(
      base_ode, y0_d, theta_v, x, x_int, &msgs);

  std::vector<double> z(z_size, 0);
  std::vector<double> dz_dt(z_size, 0);
  double t = 10;

  EXPECT_TRUE(stan::math::empty_nested());
  EXPECT_NO_THROW(coupled_system_dv(z, dz_dt, t));
  EXPECT_TRUE(stan::math::empty_nested());
}

TEST_F(StanAgradRevOde, memory_recovery_exception_dv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  std::string message = "ode throws";

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * M;

  for (size_t n = 0; n < N + 1; n++) {
    std::stringstream scoped_message;
    scoped_message << "iteration " << n;
    SCOPED_TRACE(scoped_message.str());
    mock_throwing_ode_functor<std::logic_error> throwing_ode(message, 1);

    std::vector<double> y0_d(N, 0.0);
    std::vector<var> theta_v(M, 0.0);

    coupled_ode_system<mock_throwing_ode_functor<std::logic_error>, double, var>
        coupled_system_dv(throwing_ode, y0_d, theta_v, x, x_int, &msgs);

    std::vector<double> z(z_size, 0);
    std::vector<double> dz_dt(z_size, 0);
    double t = 10;

    EXPECT_TRUE(stan::math::empty_nested());
    EXPECT_THROW_MSG(coupled_system_dv(z, dz_dt, t), std::logic_error, message);
    EXPECT_TRUE(stan::math::empty_nested());
  }
}

// ******************** VD ****************************

TEST_F(StanAgradRevOde, coupled_ode_system_vd) {
  using stan::math::coupled_ode_system;

  stan::math::start_nested();

  harm_osc_ode_fun harm_osc;

  std::vector<double> theta;
  std::vector<double> z0;
  std::vector<stan::math::var> y0_var;
  std::vector<double> y0_adj;
  double t0;
  const size_t N = 2;
  const size_t z_size = N + N * N;
  std::vector<double> dz_dt(z_size, 0);

  double gamma(0.15);
  t0 = 0;

  theta.push_back(gamma);

  z0.push_back(1.0);
  z0.push_back(0.5);
  z0.push_back(1.0);
  z0.push_back(0.0);
  z0.push_back(0.0);
  z0.push_back(1.0);

  y0_var.push_back(1.0);
  y0_var.push_back(0.5);

  std::size_t stack_size = stan::math::nested_size();

  coupled_ode_system<harm_osc_ode_fun, stan::math::var, double> system(
      harm_osc, y0_var, theta, x, x_int, &msgs);

  EXPECT_EQ(stack_size, stan::math::nested_size())
      << "expecting no new things on the stack";

  system(z0, dz_dt, t0);

  EXPECT_FLOAT_EQ(0.5, dz_dt[0]);
  EXPECT_FLOAT_EQ(-1.0 - 0.15 * 0.5, dz_dt[1]);
  EXPECT_FLOAT_EQ(0.0 * 1.0 + 1.0 * 0.0, dz_dt[2]);
  EXPECT_FLOAT_EQ(-1.0 * 1.0 - 0.15 * 0.0, dz_dt[3]);
  EXPECT_FLOAT_EQ(0.0 * 0.0 + 1.0 * 1.0, dz_dt[4]);
  EXPECT_FLOAT_EQ(-1.0 * 0.0 - 0.15 * 1.0, dz_dt[5]);

  stan::math::recover_memory_nested();
}
TEST_F(StanAgradRevOde, initial_state_vd) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;

  std::vector<var> y0_v(N, 0.0);
  std::vector<double> theta_d(M, 0.0);

  for (size_t n = 0; n < N; n++)
    y0_v[n] = n + 1;
  for (size_t m = 0; m < M; m++)
    theta_d[m] = 10 * (m + 1);

  coupled_ode_system<mock_ode_functor, var, double> coupled_system_vd(
      base_ode, y0_v, theta_d, x, x_int, &msgs);

  std::vector<double> state;

  state = coupled_system_vd.initial_state();
  for (size_t n = 0; n < N; n++)
    EXPECT_FLOAT_EQ(n + 1, state[n]);
  for (size_t i = 0; i < N; i++)
    for (size_t j = 0; j < N; j++)
      EXPECT_FLOAT_EQ(i == j ? 1.0 : 0.0, state[N + i + j * N]);
}
TEST_F(StanAgradRevOde, size_vd) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * N;

  std::vector<var> y0_v(N, 0.0);
  std::vector<double> theta_d(M, 0.0);

  coupled_ode_system<mock_ode_functor, var, double> coupled_system_vd(
      base_ode, y0_v, theta_d, x, x_int, &msgs);

  EXPECT_EQ(z_size, coupled_system_vd.size());
}

TEST_F(StanAgradRevOde, memory_recovery_vd) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * N;

  std::vector<var> y0_v(N, 0.0);
  std::vector<double> theta_d(M, 0.0);

  coupled_ode_system<mock_ode_functor, var, double> coupled_system_vd(
      base_ode, y0_v, theta_d, x, x_int, &msgs);

  std::vector<double> z(z_size, 0);
  std::vector<double> dz_dt(z_size, 0);
  double t = 10;

  EXPECT_TRUE(stan::math::empty_nested());
  EXPECT_NO_THROW(coupled_system_vd(z, dz_dt, t));
  EXPECT_TRUE(stan::math::empty_nested());
}

TEST_F(StanAgradRevOde, memory_recovery_exception_vd) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  std::string message = "ode throws";

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * N;

  for (size_t n = 0; n < N + 1; n++) {
    std::stringstream scoped_message;
    scoped_message << "iteration " << n;
    SCOPED_TRACE(scoped_message.str());
    mock_throwing_ode_functor<std::logic_error> throwing_ode(message, 1);

    std::vector<var> y0_v(N, 0.0);
    std::vector<double> theta_d(M, 0.0);

    coupled_ode_system<mock_throwing_ode_functor<std::logic_error>, var, double>
        coupled_system_vd(throwing_ode, y0_v, theta_d, x, x_int, &msgs);

    std::vector<double> z(z_size, 0);
    std::vector<double> dz_dt(z_size, 0);
    double t = 10;

    EXPECT_TRUE(stan::math::empty_nested());
    EXPECT_THROW_MSG(coupled_system_vd(z, dz_dt, t), std::logic_error, message);
    EXPECT_TRUE(stan::math::empty_nested());
  }
}

// ******************** VV ****************************

TEST_F(StanAgradRevOde, coupled_ode_system_vv) {
  using stan::math::coupled_ode_system;

  stan::math::start_nested();
  const size_t N = 2;
  const size_t M = 1;
  const size_t z_size = N + N * N + N * M;

  std::vector<stan::math::var> y0_var;
  y0_var.push_back(1.0);
  y0_var.push_back(0.5);

  std::vector<stan::math::var> theta_var;
  theta_var.push_back(0.15);

  harm_osc_ode_fun harm_osc;

  std::size_t stack_size = stan::math::nested_size();

  coupled_ode_system<harm_osc_ode_fun, stan::math::var, stan::math::var> system(
      harm_osc, y0_var, theta_var, x, x_int, &msgs);

  EXPECT_EQ(stack_size, stan::math::nested_size())
      << "expecting no new things on the stack";

  std::vector<double> z0(z_size, 0);
  z0[0] = 1.0;
  z0[1] = 0.5;
  z0[2] = 1.0;
  z0[5] = 1.0;

  double t0;
  t0 = 0;

  std::vector<double> dz_dt(z_size);
  system(z0, dz_dt, t0);

  std::vector<double> y0_double(2);
  y0_double[0] = 1.0;
  y0_double[1] = 0.5;

  std::vector<double> theta_double(1);
  theta_double[0] = 0.15;

  std::vector<double> dy_dt_base
      = harm_osc(0.0, y0_double, theta_double, x, x_int, &msgs);

  EXPECT_FLOAT_EQ(dy_dt_base[0], dz_dt[0]);
  EXPECT_FLOAT_EQ(dy_dt_base[1], dz_dt[1]);
  EXPECT_FLOAT_EQ(0, dz_dt[2]);
  EXPECT_FLOAT_EQ(-1, dz_dt[3]);
  EXPECT_FLOAT_EQ(1, dz_dt[4]);
  EXPECT_FLOAT_EQ(-0.15, dz_dt[5]);
  EXPECT_FLOAT_EQ(0, dz_dt[6]);
  EXPECT_FLOAT_EQ(-0.5, dz_dt[7]);

  stan::math::recover_memory_nested();
}
TEST_F(StanAgradRevOde, initial_state_vv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;

  std::vector<var> y0_v(N, 0.0);
  std::vector<var> theta_v(M, 0.0);

  for (size_t n = 0; n < N; n++)
    y0_v[n] = n + 1;
  for (size_t m = 0; m < M; m++)
    theta_v[m] = 10 * (m + 1);

  coupled_ode_system<mock_ode_functor, var, var> coupled_system_vv(
      base_ode, y0_v, theta_v, x, x_int, &msgs);

  std::vector<double> state = coupled_system_vv.initial_state();
  for (size_t n = 0; n < N; n++)
    EXPECT_FLOAT_EQ(n + 1, state[n]);
  for (size_t i = 0; i < N; i++)
    for (size_t j = 0; j < N; j++)
      EXPECT_FLOAT_EQ(i == j ? 1.0 : 0.0, state[N + i + j * N]);
  for (size_t n = N + N * N; n < N + N * N + N * M; n++)
    EXPECT_FLOAT_EQ(0.0, state[n]);
}
TEST_F(StanAgradRevOde, size_vv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * N + N * M;

  std::vector<var> y0_v(N, 0.0);
  std::vector<var> theta_v(M, 0.0);

  coupled_ode_system<mock_ode_functor, var, var> coupled_system_vv(
      base_ode, y0_v, theta_v, x, x_int, &msgs);

  EXPECT_EQ(z_size, coupled_system_vv.size());
}

TEST_F(StanAgradRevOde, memory_recovery_vv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  mock_ode_functor base_ode;

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * N + N * M;

  std::vector<var> y0_v(N, 0.0);
  std::vector<var> theta_v(M, 0.0);

  coupled_ode_system<mock_ode_functor, var, var> coupled_system_vv(
      base_ode, y0_v, theta_v, x, x_int, &msgs);

  std::vector<double> z(z_size, 0);
  std::vector<double> dz_dt(z_size, 0);
  double t = 10;

  EXPECT_TRUE(stan::math::empty_nested());
  EXPECT_NO_THROW(coupled_system_vv(z, dz_dt, t));
  EXPECT_TRUE(stan::math::empty_nested());
}

TEST_F(StanAgradRevOde, memory_recovery_exception_vv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;
  std::string message = "ode throws";

  const size_t N = 3;
  const size_t M = 4;
  const size_t z_size = N + N * N + N * M;

  for (size_t n = 0; n < N + 1; n++) {
    std::stringstream scoped_message;
    scoped_message << "iteration " << n;
    SCOPED_TRACE(scoped_message.str());
    mock_throwing_ode_functor<std::logic_error> throwing_ode(message, 1);

    std::vector<var> y0_v(N, 0.0);
    std::vector<var> theta_v(M, 0.0);

    coupled_ode_system<mock_throwing_ode_functor<std::logic_error>, var, var>
        coupled_system_vv(throwing_ode, y0_v, theta_v, x, x_int, &msgs);

    std::vector<double> z(z_size, 0);
    std::vector<double> dz_dt(z_size, 0);
    double t = 10;

    EXPECT_TRUE(stan::math::empty_nested());
    EXPECT_THROW_MSG(coupled_system_vv(z, dz_dt, t), std::logic_error, message);
    EXPECT_TRUE(stan::math::empty_nested());
  }
}
