#include <gtest/gtest.h>
#include <stan/math/prim/scal/fun/F32.hpp>

TEST(MathPrimScalFun, F32) {
  // Only testing this with z = 1
  EXPECT_NEAR(2.5,
              stan::math::F32(1.0, 1.0, 1.0, 1.0, 1.0, 0.6),
              1e-6);
  EXPECT_NEAR(1.0,
              stan::math::F32(1.0, 3.0, -1.0, 0.0, 1.0, 1.0),
              1e-6);
  // EXPECT_NEAR(0.4,
  //             stan::math::F32(1.0, 3.0, -1.0, 10.0, 0.5, 1.0),
  //             1e-6);
  // EXPECT_NEAR(1.56474718,
  //             stan::math::F32(1.0, 2.0, 3.0, 4.0, 5.0, 1.0),
  //             1e-6);
  // EXPECT_NEAR(0.73,
  //             stan::math::F32(1.0, 4.0, 1.0, -5.0, 4.0, 1.0),
  //             1e-6);
}
