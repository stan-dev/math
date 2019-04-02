#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathPrimMat, constraint_tolerance_exists) {
  using stan::math::CONSTRAINT_TOLERANCE;
  EXPECT_FLOAT_EQ(CONSTRAINT_TOLERANCE, 1E-8);
}
