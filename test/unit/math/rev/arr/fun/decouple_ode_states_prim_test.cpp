#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>


TEST(StanMathDecoupleOdeStates, decouple_ode_states_dd) {
  using stan::math::decouple_ode_states;

  std::vector<double> y0(2);
  y0[0] = 1.0;
  y0[1] = 0.5;

  std::vector<double> theta(1);
  theta[0] = 0.15;

  int N = y0.size();
  int T = 10;
  int k = 0;
  std::vector<std::vector<double> > ys_coupled(T);
  for (int t = 0; t < T; t++) {
    std::vector<double> coupled_state(N, 0.0);
    for (int n = 0; n < N; n++)
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

