#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(prob_transform, offset_multiplier_matrix_j) {
  double lp = -17.0;
  Eigen::VectorXd mu(2);
  mu << 2.0, 3.0;
  Eigen::VectorXd sigma(2);
  sigma << 5.0, 6.0;
  Eigen::VectorXd x(2);
  x << -1.0, 2.0;
  using stan::math::log;
  using stan::math::offset_multiplier_constrain;
  using stan::math::size;
  using stan::math::sum;
  EXPECT_MATRIX_EQ(mu.array() + sigma.array() * x.array(),
                   offset_multiplier_constrain(x, mu, sigma, lp));
  EXPECT_FLOAT_EQ(-17.0 + sum(log(sigma)), lp);
}

TEST(ProbTransform, offset_multiplier_matrix_Exception) {
  using stan::math::offset_multiplier_constrain;
  using stan::math::offset_multiplier_free;
  Eigen::VectorXd x(2);
  x << -1.0, 2.0;
  Eigen::VectorXd mu(2);
  mu << 2.0, 3.0;
  Eigen::VectorXd sigma(2);
  sigma << 5.0, 6.0;
  Eigen::VectorXd sigma_zero(2);
  sigma_zero << 5.0, 0.0;
  EXPECT_THROW(offset_multiplier_constrain(x, mu, sigma_zero),
               std::domain_error);
  Eigen::VectorXd mu_inf(2);
  mu_inf << 2.0, std::numeric_limits<double>::infinity();
  EXPECT_THROW(offset_multiplier_constrain(x, mu_inf, sigma),
               std::domain_error);

  Eigen::VectorXd mu_nan(2);
  mu_nan << 2.0, NAN;

  EXPECT_THROW(offset_multiplier_constrain(x, mu_nan, sigma),
               std::domain_error);
  double lp = 12;
  EXPECT_THROW(offset_multiplier_constrain(x, mu, sigma_zero, lp),
               std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(x, mu_inf, sigma, lp),
               std::domain_error);
  EXPECT_THROW(offset_multiplier_constrain(x, mu_nan, sigma, lp),
               std::domain_error);
}
TEST(prob_transform, offset_multiplier_matrix_f) {
  using stan::math::offset_multiplier_constrain;
  using stan::math::offset_multiplier_free;
  Eigen::VectorXd mu(2);
  mu << 2.0, 3.0;
  Eigen::VectorXd sigma(2);
  sigma << 5.0, 6.0;
  Eigen::VectorXd x(2);
  x << -1.0, 2.0;
  EXPECT_MATRIX_EQ(x, offset_multiplier_constrain(
                          offset_multiplier_free(x, mu, sigma), mu, sigma));
  EXPECT_MATRIX_EQ(
      x, offset_multiplier_free(offset_multiplier_constrain(x, mu, sigma), mu,
                                sigma));
}
TEST(prob_transform, offset_multiplier_matrix_f_exception) {
  using stan::math::offset_multiplier_constrain;
  using stan::math::offset_multiplier_free;
  Eigen::VectorXd mu(2);
  mu << 2.0, 3.0;
  Eigen::VectorXd sigma(2);
  sigma << 5.0, 6.0;
  Eigen::VectorXd x(2);
  x << -1.0, 2.0;
  Eigen::VectorXd sigma_zero(2);
  sigma_zero << 5.0, 0.0;
  EXPECT_THROW(offset_multiplier_free(x, mu, sigma_zero), std::domain_error);
  Eigen::VectorXd mu_inf(2);
  mu_inf << 2.0, std::numeric_limits<double>::infinity();
  EXPECT_THROW(offset_multiplier_free(x, mu_inf, sigma), std::domain_error);
  Eigen::VectorXd mu_nan(2);
  mu_nan << 2.0, NAN;
  EXPECT_THROW(offset_multiplier_free(x, mu_nan, sigma), std::domain_error);
}
