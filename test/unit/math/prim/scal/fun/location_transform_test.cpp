#include <stan/math/prim/scal.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, location) {
  EXPECT_FLOAT_EQ(2.0 - 5.0 * 1.0,
                  stan::math::location_constrain(-1.0, 2.0, 5.0));

  EXPECT_FLOAT_EQ(1.7, stan::math::location_constrain(1.7, 0, 1));
}
TEST(prob_transform, location_j) {
  double lp = -17.0;
  double L = 2.0;
  double U = 5.0;
  double x = -1.0;
  EXPECT_FLOAT_EQ(L + U * x, stan::math::location_constrain(x, L, U, lp));
  EXPECT_FLOAT_EQ(-17.0 + log(U), lp);

  double lp1 = -12.9;
  EXPECT_FLOAT_EQ(1.7, stan::math::location_constrain(1.7, 0, 1, lp1));
  EXPECT_FLOAT_EQ(-12.9, lp1);
}
TEST(ProbTransform, locationException) {
  using stan::math::location_constrain;
  EXPECT_THROW(location_constrain(5.0, 1.0, 0.0), std::domain_error);
  EXPECT_NO_THROW(location_constrain(5.0, 1.0, 0.01));
  double lp = 12;
  EXPECT_THROW(location_constrain(5.0, 1.0, 0.0, lp), std::domain_error);
  EXPECT_NO_THROW(location_constrain(5.0, 1.0, 0.01, lp));
}
TEST(prob_transform, location_f) {
  double L = -10.0;
  double U = 27.0;
  double y = 3.0;
  EXPECT_FLOAT_EQ(y, stan::math::location_constrain(
                         stan::math::location_free(y, L, U), L, U));
  EXPECT_FLOAT_EQ(y, stan::math::location_free(
                         stan::math::location_constrain(y, L, U), L, U));
  L = 0.0;
  U = 1.0;
  y = 3.0;
  EXPECT_FLOAT_EQ(y, stan::math::location_constrain(
                         stan::math::location_free(y, L, U), L, U));
  EXPECT_FLOAT_EQ(y, stan::math::location_free(
                         stan::math::location_constrain(y, L, U), L, U));
}
TEST(prob_transform, location_f_exception) {
  double L = -10.0;
  double U = -27.0;
  EXPECT_THROW(stan::math::location_free(L - 0.01, L, U), std::domain_error);
}
