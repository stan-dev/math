#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <test/unit/math/prim/fun/beta_reference.hpp>
#include <gtest/gtest.h>
#include <cmath>

TEST(MathFunctions, beta) {
  using stan::math::beta;

  EXPECT_FLOAT_EQ(beta(2.15, 1.71), 0.1936023967178879658641281697269);
  EXPECT_FLOAT_EQ(beta(7.62, 10.15), 0.0000065340071071564116445887286);
}

TEST(MathFunctions, beta_nan) {
  using stan::math::beta;
  using stan::math::INFTY;
  using stan::math::NOT_A_NUMBER;

  EXPECT_TRUE(std::isnan(beta(NOT_A_NUMBER, 2.16)));
  EXPECT_TRUE(std::isnan(beta(1.65, NOT_A_NUMBER)));

  EXPECT_TRUE(std::isnan(beta(INFTY, 2.16)));
  EXPECT_TRUE(std::isnan(beta(1.65, INFTY)));
}

TEST(MathFunctions, beta_large_arguments) {
  using stan::math::beta;
  for (const auto& r : stan::test::beta_large_reference) {
    SCOPED_TRACE(::testing::Message() << "a = " << r[0] << ", b = " << r[1]);
    EXPECT_NEAR(r[2], beta(r[0], r[1]), 1e-13 * r[2]);
  }
}

TEST(MathFunctions, beta_small_arguments_unchanged) {
  using stan::math::beta;
  using stan::math::lgamma;
  // Below 10, lbeta sums the same three lgamma values.
  const double xs[] = {1e-3, 0.37, 1.0, 2.5, 7.25, 9.99};
  for (double a : xs) {
    for (double b : xs) {
      EXPECT_EQ(std::exp(lgamma(a) + lgamma(b) - lgamma(a + b)), beta(a, b));
    }
  }
}

TEST(MathFunctions, beta_vec) {
  auto f
      = [](const auto& x1, const auto& x2) { return stan::math::beta(x1, x2); };

  Eigen::VectorXd in1 = Eigen::VectorXd::Random(6);
  Eigen::VectorXd in2 = Eigen::VectorXd::Random(6);

  stan::test::binary_scalar_tester(f, in1, in2);
}
