#include <stan/math/prim/functor/integrate_ode_std_vector_interface_adapter.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <vector>

TEST(StanMath, check_values) {
  harm_osc_ode_data_fun harm_osc;
  stan::math::internal::integrate_ode_std_vector_interface_adapter<
      harm_osc_ode_data_fun>
      harm_osc_adapted(harm_osc);

  std::vector<double> theta = {0.15};
  std::vector<double> y = {1.0, 0.5};

  std::vector<double> x(3, 1);
  std::vector<int> x_int(2, 0);

  double t = 1.0;

  Eigen::VectorXd out1
      = stan::math::to_vector(harm_osc(t, y, theta, x, x_int, nullptr));
  Eigen::VectorXd out2
      = harm_osc_adapted(t, stan::math::to_vector(y), nullptr, theta, x, x_int);

  EXPECT_MATRIX_FLOAT_EQ(out1, out2);
}
