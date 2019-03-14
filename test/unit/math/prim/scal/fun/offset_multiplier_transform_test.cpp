#include <stan/math/prim/scal.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, offset_multiplier) {
  EXPECT_FLOAT_EQ(2.0 - 5.0 * 1.0,
                  stan::math::offset_multiplier_constrain(-1.0, 2.0, 5.0));

  EXPECT_FLOAT_EQ(1.7, stan::math::offset_multiplier_constrain(1.7, 0, 1));
}
TEST(prob_transform, offset_multiplier_j) {
  double lp = -17.0;
  double L = 2.0;
  double U = 5.0;
  double x = -1.0;
  EXPECT_FLOAT_EQ(L + U * x,
                  stan::math::offset_multiplier_constrain(x, L, U, lp));
  EXPECT_FLOAT_EQ(-17.0 + log(U), lp);

  double lp1 = -12.9;
  EXPECT_FLOAT_EQ(1.7, stan::math::offset_multiplier_constrain(1.7, 0, 1, lp1));
  EXPECT_FLOAT_EQ(-12.9, lp1);
}
TEST(ProbTransform, offset_multiplierException) {
  using stan::math::offset_multiplier_constrain;
  using stan::math::offset_multiplier_free;
  EXPECT_THROW(offset_multiplier_constrain(5.0, 1.0, 0.0), std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(
                   5.0, std::numeric_limits<double>::infinity(), 1.0),
               std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(5.0, NAN, 1.0), std::domain_error);
  EXPECT_NO_THROW(offset_multiplier_constrain(5.0, 1.0, 0.01));
  EXPECT_THROW(offset_multiplier_free(5.0, 1.0, 0.0), std::domain_error);
  EXPECT_THROW(
      offset_multiplier_free(5.0, std::numeric_limits<double>::infinity(), 1.0),
      std::domain_error);
  EXPECT_THROW(offset_multiplier_free(5.0, NAN, 1.0), std::domain_error);
  EXPECT_NO_THROW(offset_multiplier_free(5.0, 1.0, 0.01));
  double lp = 12;
  EXPECT_THROW(offset_multiplier_constrain(5.0, 1.0, 0.0, lp),
               std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(
                   5.0, std::numeric_limits<double>::infinity(), 1.0, lp),
               std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(5.0, NAN, 1.0, lp),
               std::domain_error);
  EXPECT_NO_THROW(offset_multiplier_constrain(5.0, 1.0, 0.01, lp));
}
TEST(prob_transform, offset_multiplier_f) {
  double L = -10.0;
  double U = 27.0;
  double y = 3.0;
  EXPECT_FLOAT_EQ(y, stan::math::offset_multiplier_constrain(
                         stan::math::offset_multiplier_free(y, L, U), L, U));
  EXPECT_FLOAT_EQ(y,
                  stan::math::offset_multiplier_free(
                      stan::math::offset_multiplier_constrain(y, L, U), L, U));
  L = 0.0;
  U = 1.0;
  y = 3.0;
  EXPECT_FLOAT_EQ(y, stan::math::offset_multiplier_constrain(
                         stan::math::offset_multiplier_free(y, L, U), L, U));
  EXPECT_FLOAT_EQ(y,
                  stan::math::offset_multiplier_free(
                      stan::math::offset_multiplier_constrain(y, L, U), L, U));
}
TEST(prob_transform, offset_multiplier_f_exception) {
  double L = -10.0;
  double U = -27.0;
  EXPECT_THROW(stan::math::offset_multiplier_free(L - 0.01, L, U),
               std::domain_error);
}
