#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsPoissonGammaLCDF, values) {
  using stan::math::gamma_poisson_lcdf;

  // Reference values calculated by extraDistr::pgpois
  EXPECT_FLOAT_EQ(gamma_poisson_lcdf(5, 1, 6), -8.49989587616611e-06);
  EXPECT_FLOAT_EQ(gamma_poisson_lcdf(7, 8, 3), -0.0174512291323641);
  EXPECT_FLOAT_EQ(gamma_poisson_lcdf(1, 6, 1), -2.77258872223978);
  EXPECT_FLOAT_EQ(gamma_poisson_lcdf(4, 0.1, 0.1), -0.0694826332112239);

  std::vector<int> y{5, 7, 1, 4};

  Eigen::VectorXd alpha(4);
  alpha << 1, 8, 6, 0.1;

  Eigen::VectorXd beta(4);
  beta << 6, 3, 1, 0.1;

  double cdf = -8.49989587616611e-06 - 0.0174512291323641 - 2.77258872223978
               - 0.0694826332112239;

  EXPECT_FLOAT_EQ(gamma_poisson_lcdf(y, alpha, beta), cdf);
}
TEST(ProbDistributionsPoissonGammaLCDF, errors) {
  EXPECT_NO_THROW(stan::math::gamma_poisson_lcdf(0, 6, 5));
  EXPECT_NO_THROW(stan::math::gamma_poisson_lcdf(10, 0.5, 15));
  EXPECT_NO_THROW(stan::math::gamma_poisson_lcdf(31, 2, 1e9));

  EXPECT_THROW(stan::math::gamma_poisson_lcdf(4, 0.1, -2), std::domain_error);
  EXPECT_THROW(stan::math::gamma_poisson_lcdf(5, 6, -2), std::domain_error);
  EXPECT_THROW(stan::math::gamma_poisson_lcdf(6, -6, -0.1), std::domain_error);
  EXPECT_THROW(
      stan::math::gamma_poisson_lcdf(7, stan::math::positive_infinity(), 2),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gamma_poisson_lcdf(8, stan::math::positive_infinity(), 6),
      std::domain_error);
  EXPECT_THROW(
      stan::math::gamma_poisson_lcdf(9, 2, stan::math::positive_infinity()),
      std::domain_error);
}
