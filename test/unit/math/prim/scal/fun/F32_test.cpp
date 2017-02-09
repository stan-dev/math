#include <gtest/gtest.h>
#include <stan/math/prim/scal/fun/F32.hpp>

TEST(MathPrimScalFun, F32) {
  EXPECT_NEAR(2.5,
              stan::math::F32(1.0, 1.0, 1.0, 1.0, 1.0, 0.6),
              1e-6);
  EXPECT_NEAR(11.28915378492300834453857665243661995978358572684,
              stan::math::F32(1.0, 31.0, -27.0, 19.0, -41.0, 1.0),
              1e-6);
  EXPECT_NEAR(-0.2,
              stan::math::F32(1.0, 12.0, -1.0, 10.0, 1.0, 1.0),
              1e-6);
//  EXPECT_NEAR(0.4,
//              stan::math::F32(1.0, 3.0, -1.0, 10.0, 0.5, 1.0),
//              1e-6);
//  EXPECT_NEAR(1.56474718,
//              stan::math::F32(1.0, 2.0, 3.0, 4.0, 5.0, 1.0),
//              1e-6);
//  EXPECT_NEAR(1.0,
//              stan::math::F32(1.0, 4.0, 0.0, -5.0, 4.0, 1.0),
//              1e-6);
}


