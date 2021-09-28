#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, apply_scalar_binary_copy_scalars) {
  int n = 1;
  Eigen::VectorXd lambda(1);
  lambda << 0.4;

  const auto& y1 = stan::math::apply_scalar_binary(
      n + 1, lambda,
      [&](const auto& c, const auto& d) { return stan::math::gamma_q(c, d); });
  EXPECT_NO_THROW(stan::math::sum(y1));

  const auto& y2 = stan::math::apply_scalar_binary(
      lambda, n + 1,
      [&](const auto& d, const auto& c) { return stan::math::gamma_q(c, d); });
  EXPECT_NO_THROW(stan::math::sum(y2));
}
