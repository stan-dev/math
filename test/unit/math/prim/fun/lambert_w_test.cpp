#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

TEST(MathFunctions, lambert_w) {
  using stan::math::exp;
  using stan::math::lambert_w0;
  using stan::math::lambert_wm1;

  EXPECT_FLOAT_EQ(-1.0, lambert_w0(-1 / exp(1)));
  EXPECT_FLOAT_EQ(1.7455280027406994, lambert_w0(10.));
  EXPECT_FLOAT_EQ(1.7455280027406994, lambert_w0(10));
  EXPECT_FLOAT_EQ(-1.0, lambert_wm1(-1 / exp(1)));
  EXPECT_FLOAT_EQ(lambert_wm1(-std::numeric_limits<double>::min()),
                  -714.96865723796634);
}

TEST(MathFunctions, lambert_wn1_at_0) {
  EXPECT_TRUE(std::isinf(stan::math::lambert_wm1(0)));
}

TEST(MathFunctions, lambert_w0_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << -0.2, 0, 42, 66, 1024.512;
  stan::math::multiply(a, stan::math::lambert_w0(b));
}

TEST(MathFunctions, lambert_wn1_works_with_other_functions) {
  Eigen::VectorXd a(5);
  a << 1.1, 1.2, 1.3, 1.4, 1.5;
  Eigen::RowVectorXd b(5);
  b << -0.2, -0.15, -0.1, -0.05, -0.001;
  stan::math::multiply(a, stan::math::lambert_wm1(b));
}
