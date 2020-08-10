#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/prim/functor/mock_ode_functor.hpp>
#include <vector>
#include <string>

TEST(MathPrim, ode_store_sensitivities) {
  mock_ode_functor base_ode;

  size_t N = 5;

  int ignored1 = 0;
  double ignored2 = 0.0;
  Eigen::VectorXd ignored3;
  std::vector<std::vector<double>> ignored4;

  std::vector<double> coupled_state(N);
  for (size_t i = 0; i < coupled_state.size(); ++i)
    coupled_state[i] = i + 1;

  Eigen::VectorXd output = stan::math::ode_store_sensitivities(
      base_ode, coupled_state, Eigen::VectorXd::Ones(3).eval(), 0, 0.0, nullptr,
      ignored1, ignored2, ignored3, ignored4);

  EXPECT_EQ(output.size(), coupled_state.size());
  EXPECT_MATRIX_FLOAT_EQ(output, stan::math::to_vector(coupled_state).eval());
}
