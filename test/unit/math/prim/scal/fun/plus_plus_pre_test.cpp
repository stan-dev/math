#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, plusPlusPreDoubleVal) {
  using stan::math::plus_plus_pre;
  double x = 4;
  double y = plus_plus_pre(x);
  EXPECT_FLOAT_EQ(5, x);
  EXPECT_FLOAT_EQ(5, y);
}
TEST(MathFunctions, plusPlusPreDoubleRef) {
  using stan::math::plus_plus_pre;
  double x = 4;
  EXPECT_EQ(&x, &plus_plus_pre(x));
}

TEST(MathFunctions, plusPlusPreIntVal) {
  using stan::math::plus_plus_pre;
  int x = 4;
  int y = plus_plus_pre(x);
  EXPECT_FLOAT_EQ(5, x);
  EXPECT_FLOAT_EQ(5, y);
}
TEST(MathFunctions, plusPlusPreIntRef) {
  using stan::math::plus_plus_pre;
  int x = 4;
  EXPECT_EQ(&x, &plus_plus_pre(x));
}
