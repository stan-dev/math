#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, apply_scalar_binary_copy_scalars) {
  int n = 1;
  Eigen::VectorXd lambda(1);
  lambda << 0.4;

  const auto& y1 = stan::math::apply_scalar_binary(
      [](const auto& c, const auto& d) { return stan::math::gamma_q(c, d); },
      n + 1, lambda);
  EXPECT_NO_THROW(stan::math::sum(y1));

  const auto& y2 = stan::math::apply_scalar_binary(
      [](const auto& d, const auto& c) { return stan::math::gamma_q(c, d); },
      lambda, n + 1);
  EXPECT_NO_THROW(stan::math::sum(y2));
}
