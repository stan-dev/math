#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsPoissonGammaCDF, values) {
  using stan::math::gamma_poisson_cdf;

  // Reference values calculated by extraDistr::pgpois
  EXPECT_FLOAT_EQ(gamma_poisson_cdf(5, 1, 6), 0.999991500140248);
  EXPECT_FLOAT_EQ(gamma_poisson_cdf(7, 8, 3), 0.982700161635877);
  EXPECT_FLOAT_EQ(gamma_poisson_cdf(1, 6, 1), 0.0625);
  EXPECT_FLOAT_EQ(gamma_poisson_cdf(4, 0.1, 0.1), 0.932876334310129);

  std::vector<int> y{5, 7, 1, 4};

  Eigen::VectorXd alpha(4);
  alpha << 1, 8, 6, 0.1;

  Eigen::VectorXd beta(4);
  beta << 6, 3, 1, 0.1;

  double cdf
      = 0.999991500140248 * 0.982700161635877 * 0.0625 * 0.932876334310129;

  EXPECT_FLOAT_EQ(gamma_poisson_cdf(y, alpha, beta), cdf);
}
TEST(ProbDistributionsPoissonGammaCDF, errors) {
  EXPECT_NO_THROW(stan::math::gamma_poisson_cdf(0, 6, 5));
  EXPECT_NO_THROW(stan::math::gamma_poisson_cdf(10, 0.5, 15));
  EXPECT_NO_THROW(stan::math::gamma_poisson_cdf(31, 2, 1e9));

  EXPECT_THROW(stan::math::gamma_poisson_cdf(4, 0.1, -2), std::domain_error);
  EXPECT_THROW(stan::math::gamma_poisson_cdf(5, 6, -2), std::domain_error);
  EXPECT_THROW(stan::math::gamma_poisson_cdf(6, -6, -0.1), std::domain_error);
  EXPECT_THROW(
      stan::math::gamma_poisson_cdf(7, stan::math::positive_infinity(), 2),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gamma_poisson_cdf(8, stan::math::positive_infinity(), 6),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gamma_poisson_cdf(9, 2, stan::math::positive_infinity()),
      std::domain_error);
}
