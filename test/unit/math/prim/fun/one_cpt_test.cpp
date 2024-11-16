#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

TEST(MathFunctions, one_cpt) {
  using stan::math::one_cpt;

  Eigen::VectorXd y(2), rate(2), yt(2);
  y(0) = 745;
  y(1) = 100;
  rate(0) = 1200;
  rate(1) = 200;
  const double ka = 0.0;
  const double k10 = 50.0/80.0;	// CL/V
  yt = one_cpt(y, 0.1, ka, k10, rate);
  EXPECT_FLOAT_EQ(yt(0), 865.0);
  EXPECT_FLOAT_EQ(yt(1), 113.32912);
}
