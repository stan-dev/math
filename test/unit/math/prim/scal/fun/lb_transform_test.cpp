#include <stan/math/prim/scal.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(prob_transform, lb) {
  EXPECT_FLOAT_EQ(exp(-1.0) + 2.0, stan::math::lb_constrain(-1.0,2.0));
  EXPECT_FLOAT_EQ(7.9, 
                  stan::math::lb_constrain(7.9, -std::numeric_limits<double>::infinity()));
}
TEST(prob_transform, lb_j) {
  double lp = 15.0;
  EXPECT_FLOAT_EQ(exp(-1.0) + 2.0, stan::math::lb_constrain(-1.0,2.0,lp));
  EXPECT_FLOAT_EQ(15.0 - 1.0, lp);

  double lp2 = 8.6;
  EXPECT_FLOAT_EQ(7.9, 
                  stan::math::lb_constrain(7.9, -std::numeric_limits<double>::infinity(),
                                           lp2));
  EXPECT_FLOAT_EQ(8.6, lp2);
}
TEST(prob_transform, lb_f) {
  EXPECT_FLOAT_EQ(log(3.0 - 2.0), stan::math::lb_free(3.0,2.0));
  EXPECT_FLOAT_EQ(1.7, stan::math::lb_free(1.7, -std::numeric_limits<double>::infinity()));
}
TEST(prob_transform, lb_f_exception) {
  double lb = 2.0;
  EXPECT_THROW (stan::math::lb_free(lb - 0.01, lb), std::domain_error);
}
TEST(prob_transform, lb_rt) {
  double x = -1.0;
  double xc = stan::math::lb_constrain(x,2.0);
  double xcf = stan::math::lb_free(xc,2.0);
  EXPECT_FLOAT_EQ(x,xcf);
  double xcfc = stan::math::lb_constrain(xcf,2.0);
  EXPECT_FLOAT_EQ(xc,xcfc);
}

