#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>


TEST(StanMathRevDecoupleOdeStates, decouple_ode_states_dv) {
  using stan::math::var;

  size_t T = 10;

  std::vector<double> y0(2);
  y0[0] = 1.0;
  y0[1] = 0.5;

  std::vector<double> theta_dbl(1);
  theta_dbl[0] = 0.15;
  std::vector<var> theta(theta_dbl.begin(), theta_dbl.end());

  size_t S = 1;
  size_t N = 2;
  size_t size = N * (1+S);
  size_t k = 0;
  std::vector<std::vector<double> > ys_coupled(T);
  for (size_t t = 0; t < T; t++) {
    std::vector<double> coupled_state(size, 0.0);
    for (size_t n = 0; n < size; n++)
      coupled_state[n] = ++k;
    ys_coupled[t] = coupled_state;
  }

  std::vector<std::vector<var> > ys;
  ys = decouple_ode_states(ys_coupled, y0, theta);

  ASSERT_EQ(T, ys.size());
  for (size_t t = 0; t < T; t++)
    ASSERT_EQ(2U, ys[t].size());

  for (size_t t = 0; t < T; t++)
    for (size_t n = 0; n < 2; n++)
      EXPECT_FLOAT_EQ(ys_coupled[t][n], ys[t][n].val());
}

TEST(StanMathRevDecoupleOdeStates, decouple_ode_states_vd) {
  using stan::math::var;

  size_t T = 10;

  std::vector<double> y0_d(2);
  std::vector<double> theta(1);

  y0_d[0] = 1.0;
  y0_d[1] = 0.5;
  theta[0] = 0.15;

  std::vector<var> y0_v(y0_d.begin(), y0_d.end());

  size_t N = 2;
  size_t S = N;
  size_t size = N * (1+S);
  size_t k = 0;
  std::vector<std::vector<double> > ys_coupled(T);
  for (size_t t = 0; t < T; t++) {
    std::vector<double> coupled_state(size, 0.0);
    for (size_t n = 0; n < size; n++)
      coupled_state[n] = ++k;
    ys_coupled[t] = coupled_state;
  }

  std::vector<std::vector<var> > ys;
  ys = decouple_ode_states(ys_coupled, y0_v, theta);

  ASSERT_EQ(T, ys.size());
  for (size_t t = 0; t < T; t++)
    ASSERT_EQ(2U, ys[t].size());

  // note: decouple states operation above does not do any
  // shifting. Integrator gives plain solution.
  for (size_t t = 0; t < T; t++)
    for (size_t n = 0; n < 2; n++)
      EXPECT_FLOAT_EQ(ys_coupled[t][n], // + y0_v[n].val(),
                      ys[t][n].val());
}

TEST(StanMathRevDecoupleOdeStates, decouple_ode_states_vv) {
  using stan::math::decouple_ode_states;
  using stan::math::var;

  std::vector<double> y0_d(2);
  std::vector<double> theta_d(1);

  y0_d[0] = 1.0;
  y0_d[1] = 0.5;
  theta_d[0] = 0.15;

  std::vector<var> y0_v(y0_d.begin(), y0_d.end());
  std::vector<var> theta_v(theta_d.begin(), theta_d.end());

  size_t N = 2;
  size_t S = N+1;
  size_t size = N * (1+S);
  size_t T = 10;
  size_t k = 0;
  std::vector<std::vector<double> > ys_coupled(T);
  for (size_t t = 0; t < T; t++) {
    std::vector<double> coupled_state(size, 0.0);
    for (size_t n = 0; n < size; n++)
      coupled_state[n] = ++k;
    ys_coupled[t] = coupled_state;
  }

  std::vector<std::vector<var> > ys;
  ys = decouple_ode_states(ys_coupled, y0_v, theta_v);

  ASSERT_EQ(T, ys.size());
  for (size_t t = 0; t < T; t++)
    ASSERT_EQ(2U, ys[t].size());

  // note: decouple states operation above does not do any
  // shifting. Integrator gives plain solution.
  for (size_t t = 0; t < T; t++)
    for (size_t n = 0; n < 2; n++)
      EXPECT_FLOAT_EQ(ys_coupled[t][n], // + y0[n].val(),
                      ys[t][n].val());
}

