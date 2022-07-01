#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsPoissonGamma, values) {
  using stan::math::gamma_poisson_lpmf;

  // Reference values calculated by extraDistr::dgpois
  EXPECT_FLOAT_EQ(gamma_poisson_lpmf(10, 4, 2), -6.95199150829391);
  EXPECT_FLOAT_EQ(gamma_poisson_lpmf(1, 1, 1), -1.38629436111989);
  EXPECT_FLOAT_EQ(gamma_poisson_lpmf(0, 1, 1), -0.693147180559945);
  EXPECT_FLOAT_EQ(gamma_poisson_lpmf(1000000, 1, 1), -693147.873707126);

  std::vector<int> y{10, 1, 0, 1000000};

  Eigen::VectorXd alpha(4);
  alpha << 4, 1, 1, 1;

  Eigen::VectorXd beta(4);
  beta << 2, 1, 1, 1;

  double lp = -6.95199150829391 - 1.38629436111989 - 0.693147180559945
              - 693147.873707126;

  EXPECT_FLOAT_EQ(gamma_poisson_lpmf(y, alpha, beta), lp);

  EXPECT_NO_THROW(stan::math::gamma_poisson_lpmf(1, 6, 2));
  EXPECT_NO_THROW(stan::math::gamma_poisson_lpmf(2, 0.5, 1));
  EXPECT_NO_THROW(stan::math::gamma_poisson_lpmf(3, 1e9, 1));

  EXPECT_THROW(stan::math::gamma_poisson_lpmf(4, 0, -2), std::domain_error);
  EXPECT_THROW(stan::math::gamma_poisson_lpmf(5, 6, -2), std::domain_error);
  EXPECT_THROW(stan::math::gamma_poisson_lpmf(6, -6, -0.1), std::domain_error);
  EXPECT_THROW(
      stan::math::gamma_poisson_lpmf(7, stan::math::positive_infinity(), 2),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gamma_poisson_lpmf(8, stan::math::positive_infinity(), 6),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gamma_poisson_lpmf(9, 2, stan::math::positive_infinity()),
      std::domain_error);
}
