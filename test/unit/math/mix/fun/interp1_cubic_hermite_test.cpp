#include <stan/math/mix/fun.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <tuple>

TEST(MathPrimFun, interp1CubicHermiteSetupCoefficients) {
  Eigen::VectorXd xk(4);
  Eigen::VectorXd yk(4);
  Eigen::VectorXd dydxk(4);
  xk << 0.0, 1.0, 2.5, 4.0;
  yk << 10.0, 12.0, 9.0, 11.0;
  dydxk << 1.0, 0.5, -0.5, 1.5;

  auto W = stan::math::interp1_cubic_hermite_setup(xk, yk, dydxk);
  const auto& coef = std::get<0>(W);

  // Interval [0.0, 1.0]
  EXPECT_FLOAT_EQ(10.0, coef(0,0)); // y0
  EXPECT_FLOAT_EQ(1.0, coef(0,1)); // s0
  EXPECT_FLOAT_EQ(3.5, coef(0,2)); // a2
  EXPECT_FLOAT_EQ(-2.5, coef(0,3)); // a3

  // Interval [1.0, 2.5]
  EXPECT_FLOAT_EQ(12.0, coef(1,0));
  EXPECT_FLOAT_EQ(0.5, coef(1,1));
  EXPECT_FLOAT_EQ(-4.333333333333333, coef(1,2));
  EXPECT_FLOAT_EQ(1.7777777777777777, coef(1,3));

  // Interval [2.5, 4.0]
  EXPECT_FLOAT_EQ(9.0, coef(2,0));
  EXPECT_FLOAT_EQ(-0.5, coef(2,1));
  EXPECT_FLOAT_EQ(2.333333333333333, coef(2,2));
  EXPECT_FLOAT_EQ(-0.7407407407407406, coef(2,3));
}

TEST(MathPrimFun, interp1CubicHermiteEvaluation) {
  Eigen::VectorXd xk(3);
  Eigen::VectorXd yk(3);
  Eigen::VectorXd dydxk(3);
  xk << 0.0, 1.0, 3.0;
  yk << 2.0, 4.0, 8.0;
  dydxk << 1.0, 2.0, 1.5;
  auto W = stan::math::interp1_cubic_hermite_setup(xk, yk, dydxk);

  // Test points within intervals
  EXPECT_NEAR(2.875, stan::math::interp1_cubic_hermite_eval(0.5, W), 1e-8);
  EXPECT_NEAR(4.409, stan::math::interp1_cubic_hermite_eval(1.2, W), 1e-8);
  EXPECT_NEAR(6.125, stan::math::interp1_cubic_hermite_eval(2.0, W), 1e-8);

  // Test endpoints
  EXPECT_FLOAT_EQ(2.0, stan::math::interp1_cubic_hermite_eval(0.0, W));
  EXPECT_FLOAT_EQ(8.0, stan::math::interp1_cubic_hermite_eval(3.0, W));
}

TEST(MathPrimFun, interp1CubicHermiteEdgeCases) {
  // Tightly spaced knots
  Eigen::VectorXd xk_tight(2);
  Eigen::VectorXd yk_tight(2);
  Eigen::VectorXd dydxk_tight(2);
  xk_tight << 0.0, 1e-10;
  yk_tight << 1.0, 1.0;
  dydxk_tight << 0.0, 0.0;
  auto W_tight = stan::math::interp1_cubic_hermite_setup(xk_tight, yk_tight, dydxk_tight);
  EXPECT_FLOAT_EQ(1.0, stan::math::interp1_cubic_hermite_eval(0.5e-10, W_tight));
  EXPECT_FLOAT_EQ(1.0, stan::math::interp1_cubic_hermite_eval(0.9e-10, W_tight));

  // Large values
  Eigen::VectorXd xk_large(3);
  Eigen::VectorXd yk_large(3);
  Eigen::VectorXd dydxk_large(3);
  xk_large << 1e10, 2e10, 3e10;
  yk_large << 1e20, 2e20, 3e20;
  dydxk_large << 1e15, 2e15, 3e15;
  auto W_large = stan::math::interp1_cubic_hermite_setup(xk_large, yk_large, dydxk_large);
  EXPECT_FLOAT_EQ(-1.24985e24, stan::math::interp1_cubic_hermite_eval(1.5e10, W_large));
  EXPECT_FLOAT_EQ(-2.2497028e24, stan::math::interp1_cubic_hermite_eval(2.9e10, W_large));
}

TEST(MathPrimFun, interp1CubicHermiteFunctionality) {
  Eigen::VectorXd xk(3);
  Eigen::VectorXd yk(3);
  Eigen::VectorXd dydxk(3);
  xk << 0.0, 1.0, 3.0;
  yk << 2.0, 4.0, 8.0;
  dydxk << 1.0, 2.0, 1.5;
  auto W = stan::math::interp1_cubic_hermite_setup(xk, yk, dydxk);

  // Test exact reproduction of knot values
  for (int i = 0; i < xk.size(); ++i) {
    EXPECT_FLOAT_EQ(stan::math::interp1_cubic_hermite_eval(xk(i), W), yk(i));
  }

  // Test monotonicity in appropriate regions
  EXPECT_LT(stan::math::interp1_cubic_hermite_eval(0.5, W), 
           stan::math::interp1_cubic_hermite_eval(0.6, W));
  EXPECT_LT(stan::math::interp1_cubic_hermite_eval(1.5, W),
           stan::math::interp1_cubic_hermite_eval(1.6, W));
}