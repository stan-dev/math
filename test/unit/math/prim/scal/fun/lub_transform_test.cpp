#include <stan/math/prim/scal.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, lub) {
  EXPECT_FLOAT_EQ(2.0 + (5.0 - 2.0) * stan::math::inv_logit(-1.0),
                  stan::math::lub_constrain(-1.0, 2.0, 5.0));

  EXPECT_FLOAT_EQ(1.7, stan::math::lub_constrain(
                           1.7, -std::numeric_limits<double>::infinity(),
                           +std::numeric_limits<double>::infinity()));
  EXPECT_FLOAT_EQ(stan::math::lb_constrain(1.8, 3.0),
                  stan::math::lub_constrain(
                      1.8, 3.0, +std::numeric_limits<double>::infinity()));
  EXPECT_FLOAT_EQ(stan::math::ub_constrain(1.9, -12.5),
                  stan::math::lub_constrain(
                      1.9, -std::numeric_limits<double>::infinity(), -12.5));
}
TEST(prob_transform, lub_j) {
  double lp = -17.0;
  double L = 2.0;
  double U = 5.0;
  double x = -1.0;
  EXPECT_FLOAT_EQ(L + (U - L) * stan::math::inv_logit(x),
                  stan::math::lub_constrain(x, L, U, lp));
  EXPECT_FLOAT_EQ(-17.0 + log(U - L) + log(stan::math::inv_logit(x))
                      + log(1.0 - stan::math::inv_logit(x)),
                  lp);

  double lp1 = -12.9;
  EXPECT_FLOAT_EQ(1.7, stan::math::lub_constrain(
                           1.7, -std::numeric_limits<double>::infinity(),
                           +std::numeric_limits<double>::infinity(), lp1));
  EXPECT_FLOAT_EQ(-12.9, lp1);

  double lp2 = -19.8;
  double lp2_expected = -19.8;
  EXPECT_FLOAT_EQ(stan::math::lb_constrain(1.8, 3.0, lp2_expected),
                  stan::math::lub_constrain(
                      1.8, 3.0, +std::numeric_limits<double>::infinity(), lp2));
  EXPECT_FLOAT_EQ(lp2_expected, lp2);

  double lp3 = -422;
  double lp3_expected = -422;
  EXPECT_FLOAT_EQ(
      stan::math::ub_constrain(1.9, -12.5, lp3_expected),
      stan::math::lub_constrain(1.9, -std::numeric_limits<double>::infinity(),
                                -12.5, lp3));
  EXPECT_FLOAT_EQ(lp3_expected, lp3);
}
TEST(ProbTransform, lubException) {
  using stan::math::lub_constrain;
  EXPECT_THROW(lub_constrain(5.0, 1.0, 1.0), std::domain_error);
  EXPECT_NO_THROW(lub_constrain(5.0, 1.0, 1.01));
  double lp = 12;
  EXPECT_THROW(lub_constrain(5.0, 1.0, 1.0, lp), std::domain_error);
  EXPECT_NO_THROW(lub_constrain(5.0, 1.0, 1.01, lp));
}
TEST(prob_transform, lub_f) {
  double L = -10.0;
  double U = 27.0;
  double y = 3.0;
  EXPECT_FLOAT_EQ(stan::math::logit((y - L) / (U - L)),
                  stan::math::lub_free(y, L, U));

  EXPECT_FLOAT_EQ(
      14.2, stan::math::lub_free(14.2, -std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::infinity()));
  EXPECT_FLOAT_EQ(stan::math::ub_free(-18.3, 7.6),
                  stan::math::lub_free(
                      -18.3, -std::numeric_limits<double>::infinity(), 7.6));
  EXPECT_FLOAT_EQ(stan::math::lb_free(763.9, -3122.2),
                  stan::math::lub_free(
                      763.9, -3122.2, std::numeric_limits<double>::infinity()));
}
TEST(prob_transform, lub_f_exception) {
  double L = -10.0;
  double U = 27.0;
  EXPECT_THROW(stan::math::lub_free(L - 0.01, L, U), std::domain_error);
  EXPECT_THROW(stan::math::lub_free(U + 0.01, L, U), std::domain_error);

  EXPECT_THROW(stan::math::lub_free((L + U) / 2, U, L), std::domain_error);
}
TEST(prob_transform, lub_rt) {
  double x = -1.0;
  double xc = stan::math::lub_constrain(x, 2.0, 4.0);
  double xcf = stan::math::lub_free(xc, 2.0, 4.0);
  EXPECT_FLOAT_EQ(x, xcf);
  double xcfc = stan::math::lub_constrain(xcf, 2.0, 4.0);
  EXPECT_FLOAT_EQ(xc, xcfc);
}
