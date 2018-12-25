#include <stan/math/prim/mat.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, offset_multiplier) {
  Eigen::Matrix<double, -1, 1> x(2, 1);
  x << 3, 2;
  Eigen::Matrix<double, -1, 1> mu(2, 1);
  mu << 2, 1;
  Eigen::Matrix<double, -1, -1> sigma(2, 2);
  sigma << 3, 0, 2, 1;
  Eigen::Matrix<double, -1, 1> result
      = stan::math::offset_multiplier_constrain(x, mu, sigma);
  Eigen::Matrix<double, -1, 1> expected = mu + sigma * x;
  for (size_t n = 0; n < 2; ++n) {
    EXPECT_FLOAT_EQ(result(n), expected(n));
  }

  Eigen::Matrix<double, -1, -1> sigma2(2, 2);
  sigma2 << 1, 0, 0, 1;
  result = stan::math::offset_multiplier_constrain(x, mu, sigma2);
  expected = mu + x;
  for (size_t n = 0; n < 2; ++n) {
    EXPECT_FLOAT_EQ(result(n), expected(n));
  }

  Eigen::Matrix<double, -1, 1> mu2(2, 1);
  mu2 << 0, 0;
  result = stan::math::offset_multiplier_constrain(x, mu2, sigma2);
  expected = x;
  for (size_t n = 0; n < 2; ++n) {
    EXPECT_FLOAT_EQ(result(n), expected(n));
  }
}
TEST(prob_transform, offset_multiplier_j) {
  double lp = -17.0;
  Eigen::Matrix<double, -1, 1> x(2, 1);
  x << 4, 3;
  Eigen::Matrix<double, -1, 1> mu(2, 1);
  mu << 3, 2;
  Eigen::Matrix<double, -1, -1> sigma(2, 2);
  sigma << 4, 0, 3, 2;
  Eigen::Matrix<double, -1, 1> result
      = stan::math::offset_multiplier_constrain(x, mu, sigma, lp);
  Eigen::Matrix<double, -1, 1> expected = mu + sigma * x;
  for (size_t n = 0; n < 2; ++n) {
    EXPECT_FLOAT_EQ(result(n), expected(n));
  }
  EXPECT_FLOAT_EQ(-17.0 + log(4) + log(2), lp);

  lp = -12.0;
  Eigen::Matrix<double, -1, -1> sigma2(2, 2);
  sigma2 << 1, 0, 0, 1;
  result = stan::math::offset_multiplier_constrain(x, mu, sigma2, lp);
  expected = mu + x;
  for (size_t n = 0; n < 2; ++n) {
    EXPECT_FLOAT_EQ(result(n), expected(n));
  }
  EXPECT_FLOAT_EQ(-12.0, lp);

  lp = -14.0;
  Eigen::Matrix<double, -1, 1> mu2(2, 1);
  mu2 << 0, 0;
  result = stan::math::offset_multiplier_constrain(x, mu2, sigma2, lp);
  expected = x;
  for (size_t n = 0; n < 2; ++n) {
    EXPECT_FLOAT_EQ(result(n), expected(n));
  }
  EXPECT_FLOAT_EQ(-14.0, lp);
}
/*
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
*/
