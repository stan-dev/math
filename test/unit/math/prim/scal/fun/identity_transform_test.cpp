#include <stan/math/prim/scal.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(prob_transform,identity) {
  EXPECT_FLOAT_EQ(4.0, stan::math::identity_constrain(4.0));
}
TEST(prob_transform,identity_j) {
  double lp = 1.0;
  EXPECT_FLOAT_EQ(4.0, stan::math::identity_constrain(4.0,lp));
  EXPECT_FLOAT_EQ(1.0,lp);
}
TEST(prob_transform,identity_free) {
  EXPECT_FLOAT_EQ(4.0, stan::math::identity_free(4.0));
}
TEST(prob_transform,identity_rt) {
  double x = 1.2;
  double xc = stan::math::identity_constrain(x);
  double xcf = stan::math::identity_free(xc);
  EXPECT_FLOAT_EQ(x,xcf);

  double y = -1.0;
  double yf = stan::math::identity_free(y);
  double yfc = stan::math::identity_constrain(yf);
  EXPECT_FLOAT_EQ(y,yfc);
}
