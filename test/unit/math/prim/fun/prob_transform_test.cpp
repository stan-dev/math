#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, prob) {
  EXPECT_FLOAT_EQ(stan::math::inv_logit(-1.0),
                  stan::math::prob_constrain(-1.0));
}
TEST(prob_transform, prob_j) {
  double lp = -17.0;
  double L = 0.0;
  double U = 1.0;
  double x = -1.0;
  EXPECT_FLOAT_EQ(L + (U - L) * stan::math::inv_logit(x),
                  stan::math::prob_constrain(x, lp));
  EXPECT_FLOAT_EQ(-17.0 + log(U - L) + log(stan::math::inv_logit(x))
                      + log(1.0 - stan::math::inv_logit(x)),
                  lp);
}
TEST(prob_transform, prob_f) {
  double L = 0.0;
  double U = 1.0;
  double y = 0.4;
  EXPECT_FLOAT_EQ(stan::math::logit((y - L) / (U - L)),
                  stan::math::prob_free(y));
}
TEST(prob_transform, prob_f_exception) {
  EXPECT_THROW(stan::math::prob_free(1.1), std::domain_error);
  EXPECT_THROW(stan::math::prob_free(-0.1), std::domain_error);
}
TEST(prob_transform, prob_rt) {
  double x = -1.0;
  double xc = stan::math::prob_constrain(x);
  double xcf = stan::math::prob_free(xc);
  EXPECT_FLOAT_EQ(x, xcf);
  double xcfc = stan::math::prob_constrain(xcf);
  EXPECT_FLOAT_EQ(xc, xcfc);
}
