#include <stan/math/prim/scal.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(prob_transform, corr) {
  EXPECT_FLOAT_EQ(std::tanh(-1.0), stan::math::corr_constrain(-1.0));
}
TEST(prob_transform, corr_j) {
  double lp = -17.0;
  double x = -1.0;
  EXPECT_FLOAT_EQ(std::tanh(x), stan::math::corr_constrain(x, lp));
  EXPECT_FLOAT_EQ(-17.0 + (log(1.0 - std::tanh(x) * std::tanh(x))), lp);
}
TEST(prob_transform, corr_f) {
  EXPECT_FLOAT_EQ(atanh(-0.4), 0.5 * std::log((1.0 + -0.4) / (1.0 - -0.4)));
  double y = -0.4;
  EXPECT_FLOAT_EQ(atanh(y), stan::math::corr_free(y));
}
TEST(prob_transform, corr_rt) {
  double x = -1.0;
  double xc = stan::math::corr_constrain(x);
  double xcf = stan::math::corr_free(xc);
  EXPECT_FLOAT_EQ(x, xcf);
  double xcfc = stan::math::corr_constrain(xcf);
  EXPECT_FLOAT_EQ(xc, xcfc);
}
