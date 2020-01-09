#include <gtest/gtest.h>
#include <stan/math/prim.hpp>
#include <cmath>
#include <limits>

TEST(MathPrimScalFun, FiniteDiffStepsize) {
  using stan::math::finite_diff_stepsize;
  double eps = std::numeric_limits<double>::epsilon();
  double cbrt_eps = std::cbrt(eps);
  EXPECT_FLOAT_EQ(cbrt_eps, finite_diff_stepsize(1e-20));
  EXPECT_FLOAT_EQ(cbrt_eps, finite_diff_stepsize(1e-1));
  EXPECT_FLOAT_EQ(cbrt_eps, finite_diff_stepsize(0));
  EXPECT_FLOAT_EQ(cbrt_eps * 1e10, finite_diff_stepsize(1e10));
  EXPECT_FLOAT_EQ(cbrt_eps * 1e10, finite_diff_stepsize(-1e10));
  EXPECT_FLOAT_EQ(cbrt_eps * 1e20, finite_diff_stepsize(1e20));
  EXPECT_FLOAT_EQ(cbrt_eps * 1e20, finite_diff_stepsize(-1e20));
}
