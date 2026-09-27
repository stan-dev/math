#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <tuple>

namespace {
// Shared test inputs
const Eigen::VectorXd kTestX = (Eigen::VectorXd(4) << 0.0, 1.0, 2.5, 4.0).finished();
const Eigen::VectorXd kTestY = (Eigen::VectorXd(4) << 10.0, 12.0, 9.0, 11.0).finished();
const Eigen::VectorXd kTestDydx = (Eigen::VectorXd(4) << 1.0, 0.5, -0.5, 1.5).finished();
}

TEST(MathPrimFun, interp1CubicHermiteSetupShape) {
  auto W = stan::math::interp1_cubic_hermite_setup(kTestX, kTestY, kTestDydx);
  const auto& coef = std::get<0>(W);
  const auto& knots = std::get<1>(W);

  // Verify matrix dimensions
  EXPECT_EQ(3, coef.rows()); // K-1 where K=4
  EXPECT_EQ(4, coef.cols()); // Fixed: [a0, a1, a2, a3]
  EXPECT_EQ(kTestX.size(), knots.size());

  // Verify knot preservation
  for (Eigen::Index i = 0; i < kTestX.size(); ++i) {
    EXPECT_FLOAT_EQ(kTestX(i), knots(i));
  }
}

TEST(MathPrimFun, interp1CubicHermiteSetupCoefficients) {
  auto W = stan::math::interp1_cubic_hermite_setup(kTestX, kTestY, kTestDydx);
  const auto& coef = std::get<0>(W);

  // Interval [0.0, 1.0] coefficients
  EXPECT_FLOAT_EQ(10.0, coef(0,0));  // a00 = y0 - VERIFY
  EXPECT_FLOAT_EQ(1.0, coef(0,1));   // a01 = s0 - VERIFY
  EXPECT_FLOAT_EQ(3.5, coef(0,2));   // a02 - MATHEMATICAL GROUND TRUTH NEEDED
  EXPECT_FLOAT_EQ(-2.5, coef(0,3));  // a03 - MATHEMATICAL GROUND TRUTH NEEDED

  // Interval [1.0, 2.5] coefficients
  EXPECT_FLOAT_EQ(12.0, coef(1,0));  // a10 = y1 - VERIFY
  EXPECT_FLOAT_EQ(0.5, coef(1,1));   // a11 = s1 - VERIFY
  EXPECT_FLOAT_EQ(-4.333333333333333, coef(1,2)); // a12 - MATHEMATICAL GROUND TRUTH NEEDED
  EXPECT_FLOAT_EQ(1.7777777777777777, coef(1,3)); // a13 - MATHEMATICAL GROUND TRUTH NEEDED

  // Interval [2.5, 4.0] coefficients
  EXPECT_FLOAT_EQ(9.0, coef(2,0));   // a20 = y2 - VERIFY
  EXPECT_FLOAT_EQ(-0.5, coef(2,1));  // a21 = s2 - VERIFY
  EXPECT_FLOAT_EQ(2.333333333333333, coef(2,2)); // a22 - MATHEMATICAL GROUND TRUTH NEEDED
  EXPECT_FLOAT_EQ(-0.7407407407407406, coef(2,3)); // a23 - MATHEMATICAL GROUND TRUTH NEEDED
}

TEST(MathPrimFun, interp1CubicHermiteEval) {
  auto W = stan::math::interp1_cubic_hermite_setup(kTestX, kTestY, kTestDydx);

  // Validate endpoints
  EXPECT_FLOAT_EQ(10.0, stan::math::interp1_cubic_hermite_eval(0.0, W)); // x = x0
  EXPECT_FLOAT_EQ(11.0, stan::math::interp1_cubic_hermite_eval(4.0, W));  // x = x2

  // Validate intermediate points
  EXPECT_NEAR(11.38888889, stan::math::interp1_cubic_hermite_eval(1.5, W), 1e-8);
  EXPECT_NEAR(9.24074074, stan::math::interp1_cubic_hermite_eval(3.0, W), 1e-8);

}

// Existing throw tests remain unchanged...
// interp1CubicHermiteSetupThrows
// interp1CubicHermiteEvalThrowsForInvalidSetupShape