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

TEST(ProbTransform, offset_multiplierException) {
  using stan::math::offset_multiplier_constrain;
  using stan::math::offset_multiplier_free;
  Eigen::Matrix<double, -1, 1> x(2, 1);
  x << 3, 2;
  Eigen::Matrix<double, -1, 1> mu(2, 1);
  mu << 2, 1;
  Eigen::Matrix<double, -1, -1> sigma(2, 2);
  sigma << 3, 0, 2, 1;
  Eigen::Matrix<double, -1, 1> x_bad(3, 1);
  x_bad << 3, 2, 1;
  Eigen::Matrix<double, -1, 1> x_bad2(2, 1);
  x_bad2 << NAN, 1;
  Eigen::Matrix<double, -1, 1> mu_bad(3, 1);
  mu_bad << 2, 1, 0;
  Eigen::Matrix<double, -1, 1> mu_bad2(2, 1);
  mu_bad2 << NAN, 1;
  Eigen::Matrix<double, -1, -1> sigma_bad(2, 2);
  sigma_bad << 3, -1, 2, 1;
  Eigen::Matrix<double, -1, -1> sigma_bad2(2, 2);
  sigma_bad2 << -3, 0, 2, 1;
  Eigen::Matrix<double, -1, -1> sigma_bad3(2, 2);
  sigma_bad3 << NAN, 0, 2, 1;
  
  double lp = -18.0;
 
  EXPECT_THROW(offset_multiplier_constrain(x, mu, sigma_bad), std::domain_error);

  EXPECT_THROW(offset_multiplier_constrain(x, mu, sigma_bad2), std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(x, mu, sigma_bad3), std::domain_error);

  EXPECT_THROW(offset_multiplier_constrain(x_bad, mu, sigma), std::invalid_argument);
  EXPECT_THROW(offset_multiplier_constrain(x_bad2, mu, sigma), std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(x, mu_bad, sigma), std::invalid_argument);
  EXPECT_THROW(offset_multiplier_constrain(x, mu_bad2, sigma), std::domain_error);


  EXPECT_THROW(offset_multiplier_constrain(x, mu, sigma_bad, lp), std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(x, mu, sigma_bad2, lp), std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(x, mu, sigma_bad3, lp), std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(x_bad, mu, sigma, lp), std::invalid_argument);
  EXPECT_THROW(offset_multiplier_constrain(x_bad2, mu, sigma, lp), std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(x, mu_bad, sigma, lp), std::invalid_argument);
  EXPECT_THROW(offset_multiplier_constrain(x, mu_bad2, sigma, lp), std::domain_error);

  
  EXPECT_THROW(offset_multiplier_free(x, mu, sigma_bad), std::domain_error);
  EXPECT_THROW(offset_multiplier_free(x, mu, sigma_bad2), std::domain_error);
  EXPECT_THROW(offset_multiplier_free(x, mu, sigma_bad3), std::domain_error);
  EXPECT_THROW(offset_multiplier_free(x_bad, mu, sigma), std::invalid_argument);
  EXPECT_THROW(offset_multiplier_free(x_bad2, mu, sigma), std::domain_error);
  EXPECT_THROW(offset_multiplier_free(x, mu_bad, sigma), std::invalid_argument);
  EXPECT_THROW(offset_multiplier_free(x, mu_bad2, sigma), std::domain_error);
  
}

TEST(prob_transform, offset_multiplier_f) {
  Eigen::Matrix<double, -1, 1> x(2, 1);
  x << 3, 2;
  Eigen::Matrix<double, -1, 1> mu(2, 1);
  mu << 2, 1;
  Eigen::Matrix<double, -1, -1> sigma(2, 2);
  sigma << 3, 0, 2, 1;
  Eigen::Matrix<double, -1, 1> result
      = stan::math::offset_multiplier_free(
          stan::math::offset_multiplier_constrain(x, mu, sigma), mu, sigma);
  Eigen::Matrix<double, -1, 1> expected = x;
  for (size_t n = 0; n < 2; ++n) {
    EXPECT_FLOAT_EQ(result(n), expected(n));
  }
  result
      = stan::math::offset_multiplier_constrain(
          stan::math::offset_multiplier_free(x, mu, sigma), mu, sigma);
  expected = x;
  for (size_t n = 0; n < 2; ++n) {
    EXPECT_FLOAT_EQ(result(n), expected(n));
  }
}

