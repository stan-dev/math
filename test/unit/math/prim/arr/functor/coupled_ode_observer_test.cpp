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

TEST_F(StanMathOde, observe_states_dddd) {
  using stan::math::coupled_ode_system;

  harm_osc_ode_fun harm_osc;

  std::vector<double> y0(2);
  std::vector<double> theta(1);

  y0[0] = 1.0;
  y0[1] = 0.5;
  theta[0] = 0.15;

  stan::math::coupled_ode_system<harm_osc_ode_fun, double, double>
      coupled_system(harm_osc, y0, theta, x, x_int, &msgs);

  std::vector<std::vector<double>> y;
  double t0 = 0;
  int T = 10;
  std::vector<double> ts(T);
  for (int t = 0; t < T; t++)
    ts[t] = t;

  stan::math::coupled_ode_observer<harm_osc_ode_fun, double, double, double,
                                   double>
      observer(harm_osc, y0, theta, t0, ts, x, x_int, &msgs, y);

  int k = 0;
  std::vector<std::vector<double>> ys_coupled(T);
  for (int t = 0; t < T; t++) {
    std::vector<double> coupled_state(coupled_system.size(), 0.0);
    for (int n = 0; n < coupled_system.size(); n++)
      coupled_state[n] = ++k;
    ys_coupled[t] = coupled_state;
    observer(coupled_state, ts[t]);
  }

  EXPECT_EQ(T, y.size());
  for (int t = 0; t < T; t++)
    EXPECT_EQ(2, y[t].size());

  for (int t = 0; t < T; t++)
    for (int n = 0; n < 2; n++)
      EXPECT_FLOAT_EQ(ys_coupled[t][n], y[t][n])
          << "(" << n << ", " << t << "): "
          << "for (double, double) the coupled system is the base system";

  // calling it once too often will fire an exception as the observer
  // runs out of time states.
  // note: the EXPECT_THROW_MSG calls the expression given twice which
  // is why the time-state number is 11 by the time the comparison is
  // done with the expected string (I would have expected 10 instead
  // of 11).
  std::string message
      = "coupled_ode_observer: time-state number is 10, but must be less than "
        "10";
  std::vector<double> ys_coupled_0(coupled_system.size(), 0.0);
  EXPECT_THROW_MSG(observer(ys_coupled_0, 20.), std::logic_error, message);
}
