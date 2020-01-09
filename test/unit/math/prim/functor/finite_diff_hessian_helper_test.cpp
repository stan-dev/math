#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>
#include <cmath>
#include <limits>

TEST(MathPrimMatFunctor, FiniteDiffGradientFast) {
  using stan::math::finite_diff_hessian_helper;
  auto f = [](const auto& u) { return u(0) * u(1); };
  Eigen::VectorXd x(2);
  x << -1, 2;
  EXPECT_FLOAT_EQ(0.024, finite_diff_hessian_helper(f, x, 0, 1e-3));
  EXPECT_FLOAT_EQ(-0.012, finite_diff_hessian_helper(f, x, 1, 1e-3));
}
