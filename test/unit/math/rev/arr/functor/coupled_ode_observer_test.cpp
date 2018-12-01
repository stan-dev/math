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

TEST_F(StanAgradRevOde, observe_states_dv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;

  mock_ode_functor mock_ode;

  std::vector<double> y0(2);
  y0[0] = 1.0;
  y0[1] = 0.5;

  std::vector<var> theta(1);
  theta[0] = 0.15;

  coupled_ode_system<mock_ode_functor, double, var> coupled_system(
      mock_ode, y0, theta, x, x_int, &msgs);

  std::vector<std::vector<var>> y;
  double t0 = 0;
  int T = 10;
  std::vector<double> ts(T);
  for (int t = 0; t < T; t++)
    ts[t] = t;

  stan::math::coupled_ode_observer<mock_ode_functor, double, var, double,
                                   double>
      observer(mock_ode, y0, theta, t0, ts, x, x_int, &msgs, y);

  observer(std::vector<double>(coupled_system.size(), 0.0), 0);

  size_t k = 0;
  std::vector<std::vector<double>> ys_coupled(T);
  for (size_t t = 0; t < T; t++) {
    std::vector<double> coupled_state(coupled_system.size(), 0.0);
    for (size_t n = 0; n < coupled_system.size(); n++)
      coupled_state[n] = ++k;
    ys_coupled[t] = coupled_state;
    observer(coupled_state, ts[0]);
  }

  ASSERT_EQ(T, y.size());
  for (size_t t = 0; t < T; t++)
    ASSERT_EQ(2U, y[t].size());

  for (size_t t = 0; t < T; t++) {
    for (size_t n = 0; n < 2; n++)
      EXPECT_FLOAT_EQ(ys_coupled[t][n], y[t][n].val());
    for (size_t n = 0; n < 2; n++) {
      y[t][n].grad();
      EXPECT_FLOAT_EQ(ys_coupled[t][2 + n], theta[0].adj());
      stan::math::set_zero_all_adjoints();
    }
  }
}

TEST_F(StanAgradRevOde, observe_states_vd) {
  using stan::math::coupled_ode_system;
  using stan::math::var;

  mock_ode_functor mock_ode;

  std::vector<var> y0(2);
  std::vector<double> theta(1);

  y0[0] = 1.0;
  y0[1] = 0.5;
  theta[0] = 0.15;

  coupled_ode_system<mock_ode_functor, var, double> coupled_system(
      mock_ode, y0, theta, x, x_int, &msgs);

  std::vector<std::vector<var>> y;
  double t0 = 0;
  int T = 10;
  std::vector<double> ts(T);
  for (int t = 0; t < T; t++)
    ts[t] = t;

  stan::math::coupled_ode_observer<mock_ode_functor, var, double, double,
                                   double>
      observer(mock_ode, y0, theta, t0, ts, x, x_int, &msgs, y);

  observer(std::vector<double>(coupled_system.size(), 0.0), 0);

  size_t k = 0;
  std::vector<std::vector<double>> ys_coupled(T);
  for (size_t t = 0; t < T; t++) {
    std::vector<double> coupled_state(coupled_system.size(), 0.0);
    for (size_t n = 0; n < coupled_system.size(); n++)
      coupled_state[n] = ++k;
    ys_coupled[t] = coupled_state;
    observer(coupled_state, ts[t]);
  }

  ASSERT_EQ(T, y.size());
  for (size_t t = 0; t < T; t++)
    ASSERT_EQ(2U, y[t].size());

  for (size_t t = 0; t < T; t++) {
    for (size_t n = 0; n < 2; n++) {
      EXPECT_FLOAT_EQ(ys_coupled[t][n], y[t][n].val());
    }
    for (size_t n = 0; n < 2; n++) {
      y[t][n].grad();
      EXPECT_FLOAT_EQ(ys_coupled[t][2 + n], y0[0].adj());
      EXPECT_FLOAT_EQ(ys_coupled[t][2 + 2 + n], y0[1].adj());
      stan::math::set_zero_all_adjoints();
    }
  }
}

TEST_F(StanAgradRevOde, observe_states_vv) {
  using stan::math::coupled_ode_system;
  using stan::math::var;

  harm_osc_ode_fun harm_osc;

  std::vector<var> y0(2);
  std::vector<var> theta(1);

  y0[0] = 1.0;
  y0[1] = 0.5;
  theta[0] = 0.15;

  coupled_ode_system<harm_osc_ode_fun, var, var> coupled_system(
      harm_osc, y0, theta, x, x_int, &msgs);

  std::vector<std::vector<var>> y;
  double t0 = 0;
  int T = 10;
  std::vector<double> ts(T);
  for (int t = 0; t < T; t++)
    ts[t] = t;

  stan::math::coupled_ode_observer<harm_osc_ode_fun, var, var, double, double>
      observer(harm_osc, y0, theta, t0, ts, x, x_int, &msgs, y);

  observer(std::vector<double>(coupled_system.size(), 0.0), 0);

  size_t k = 0;
  std::vector<std::vector<double>> ys_coupled(T);
  for (size_t t = 0; t < T; t++) {
    std::vector<double> coupled_state(coupled_system.size(), 0.0);
    for (size_t n = 0; n < coupled_system.size(); n++)
      coupled_state[n] = ++k;
    ys_coupled[t] = coupled_state;
    observer(coupled_state, ts[t]);
  }

  ASSERT_EQ(T, y.size());
  for (size_t t = 0; t < T; t++)
    ASSERT_EQ(2U, y[t].size());

  for (size_t t = 0; t < T; t++) {
    for (size_t n = 0; n < 2; n++)
      EXPECT_FLOAT_EQ(ys_coupled[t][n], y[t][n].val());

    for (size_t n = 0; n < 2; n++) {
      y[t][n].grad();
      EXPECT_FLOAT_EQ(ys_coupled[t][2 + n], y0[0].adj());
      EXPECT_FLOAT_EQ(ys_coupled[t][2 + 2 + n], y0[1].adj());
      EXPECT_FLOAT_EQ(ys_coupled[t][2 + 2 * 2 + n], theta[0].adj());
      stan::math::set_zero_all_adjoints();
    }
  }
}
