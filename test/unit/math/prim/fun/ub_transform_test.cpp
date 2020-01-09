#include <stan/math/prim/scal.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, ub) {
  EXPECT_FLOAT_EQ(2.0 - exp(-1.0), stan::math::ub_constrain(-1.0, 2.0));
  EXPECT_FLOAT_EQ(1.7, stan::math::ub_constrain(
                           1.7, std::numeric_limits<double>::infinity()));
}
TEST(prob_transform, ub_j) {
  double lp = 15.0;
  EXPECT_FLOAT_EQ(2.0 - exp(-1.0), stan::math::ub_constrain(-1.0, 2.0, lp));
  EXPECT_FLOAT_EQ(15.0 - 1.0, lp);

  double lp2 = 1.87;
  EXPECT_FLOAT_EQ(-5.2,
                  stan::math::ub_constrain(
                      -5.2, std::numeric_limits<double>::infinity(), lp2));
  EXPECT_FLOAT_EQ(1.87, lp2);
}
TEST(prob_transform, ub_f) {
  double y = 2.0;
  double U = 4.0;
  EXPECT_FLOAT_EQ(log(-(y - U)), stan::math::ub_free(2.0, 4.0));

  EXPECT_FLOAT_EQ(19.765, stan::math::ub_free(
                              19.765, std::numeric_limits<double>::infinity()));
}
TEST(prob_transform, ub_f_exception) {
  double ub = 4.0;
  EXPECT_THROW(stan::math::ub_free(ub + 0.01, ub), std::domain_error);
}
TEST(prob_transform, ub_rt) {
  double x = -1.0;
  double xc = stan::math::ub_constrain(x, 2.0);
  double xcf = stan::math::ub_free(xc, 2.0);
  EXPECT_FLOAT_EQ(x, xcf);
  double xcfc = stan::math::ub_constrain(xcf, 2.0);
  EXPECT_FLOAT_EQ(xc, xcfc);
}
