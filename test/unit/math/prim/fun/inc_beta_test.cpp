#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/ternary_scalar_tester.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, inc_beta) {
  using stan::math::inc_beta;

  EXPECT_FLOAT_EQ(0.0, inc_beta(0.5, 0.5, 0.0))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.333333333, inc_beta(0.5, 0.5, 0.25))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.5, inc_beta(0.5, 0.5, 0.5))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.666666667, inc_beta(0.5, 0.5, 0.75))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(1.0, inc_beta(0.5, 0.5, 1.0))
      << "reasonable values for a, b, x";

  EXPECT_FLOAT_EQ(0.0, inc_beta(0.1, 1.5, 0.0))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.9117332, inc_beta(0.1, 1.5, 0.25))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.9645342, inc_beta(0.1, 1.5, 0.5))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.9897264, inc_beta(0.1, 1.5, 0.75))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(1.0, inc_beta(0.1, 1.5, 1.0))
      << "reasonable values for a, b, x";

  EXPECT_FLOAT_EQ(0.3121373, inc_beta(0.6, 0.3, 0.5))
      << "reasonable values for a, b, x";
  EXPECT_FLOAT_EQ(0.0272, inc_beta(3, 2, 0.2))
      << "reasonable values for a, b, x";
}

TEST(MathFunctions, inc_beta_a_boundary) {
  double b = 0.5;
  double x = 0.5;

  const double inf = std::numeric_limits<double>::infinity();

  EXPECT_NO_THROW(stan::math::inc_beta(0.0, b, x));
  EXPECT_NO_THROW(stan::math::inc_beta(inf, b, x));
  EXPECT_THROW(stan::math::inc_beta(-0.01, b, x), std::domain_error);
}

TEST(MathFunctions, inc_beta_b_boundary) {
  double a = 0.5;
  double x = 0.5;

  const double inf = std::numeric_limits<double>::infinity();

  EXPECT_NO_THROW(stan::math::inc_beta(a, 0.0, x));
  EXPECT_NO_THROW(stan::math::inc_beta(inf, 1.0, x));
  EXPECT_THROW(stan::math::inc_beta(a, -0.01, x), std::domain_error);
}

TEST(MathFunctions, inc_beta_x_boundary) {
  double a = 0.5;
  double b = 0.5;

  EXPECT_NO_THROW(stan::math::inc_beta(a, b, 0.0));
  EXPECT_NO_THROW(stan::math::inc_beta(a, b, 1.0));
  EXPECT_THROW(stan::math::inc_beta(a, b, -0.01), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta(a, b, 1.01), std::domain_error);
}

TEST(MathFunctions, inc_beta_a_b_boundary) {
  double x = 0.5;

  EXPECT_THROW(stan::math::inc_beta(0.0, 0.0, x), std::domain_error);
}

TEST(MathFunctions, inc_beta_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(stan::math::inc_beta(0.5, 0.0, nan), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta(0.5, nan, 0.0), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta(0.5, nan, nan), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta(nan, 0.0, 0.0), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta(nan, 0.0, nan), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta(nan, nan, 0.0), std::domain_error);
  EXPECT_THROW(stan::math::inc_beta(nan, nan, nan), std::domain_error);
}

TEST(MathFunctions, inc_beta_vec) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    return stan::math::inc_beta(x1, x2, x3);
  };

  Eigen::VectorXd in1 = Eigen::VectorXd::Random(6).array().abs().matrix();
  Eigen::VectorXd in2 = Eigen::VectorXd::Random(6).array().abs().matrix();
  Eigen::VectorXd in3(6);
  in3 << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;

  stan::test::ternary_scalar_tester(f, in1, in2, in3);
}
