#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/arr/functor/harmonic_oscillator.hpp>
#include <test/unit/math/prim/arr/functor/mock_ode_functor.hpp>
#include <test/unit/math/prim/arr/functor/mock_throwing_ode_functor.hpp>
#include <vector>
#include <string>

struct StanMathOde : public ::testing::Test {
  std::stringstream msgs;
  std::vector<double> x;
  std::vector<int> x_int;
};

TEST_F(StanMathOde, decouple_states_dd) {
  using stan::math::coupled_ode_system;

  harm_osc_ode_fun harm_osc;

  std::vector<double> y0(2);
  std::vector<double> theta(1);

  y0[0] = 1.0;
  y0[1] = 0.5;
  theta[0] = 0.15;

  coupled_ode_system<harm_osc_ode_fun, double, double> coupled_system(
      harm_osc, y0, theta, x, x_int, &msgs);

  int T = 10;
  int k = 0;
  std::vector<std::vector<double> > ys_coupled(T);
  for (int t = 0; t < T; t++) {
    std::vector<double> coupled_state(coupled_system.size(), 0.0);
    for (int n = 0; n < coupled_system.size(); n++)
      coupled_state[n] = ++k;
    ys_coupled[t] = coupled_state;
  }

  std::vector<std::vector<double> > ys;
  ys = coupled_system.decouple_states(ys_coupled);

  ASSERT_EQ(T, ys.size());
  for (int t = 0; t < T; t++)
    ASSERT_EQ(2, ys[t].size());

  for (int t = 0; t < T; t++)
    for (int n = 0; n < 2; n++)
      EXPECT_FLOAT_EQ(ys_coupled[t][n], ys[t][n])
          << "(" << n << ", " << t << "): "
          << "for (double, double) the coupled system is the base system";
}
