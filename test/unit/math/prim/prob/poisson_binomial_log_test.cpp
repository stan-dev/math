#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(ProbDistributionsPoissonBinomial, log_matches_lpmf) {
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  int y = 2;
  vec theta(3, 1);
  theta << 0.5, 0.2, 0.7;

  EXPECT_FLOAT_EQ((stan::math::poisson_binomial_lpmf(y, theta)),
                  (stan::math::poisson_binomial_log(y, theta)));

  EXPECT_FLOAT_EQ((stan::math::poisson_binomial_lpmf<true>(y, theta)),
                  (stan::math::poisson_binomial_log<true>(y, theta)));

  EXPECT_FLOAT_EQ((stan::math::poisson_binomial_lpmf<false>(y, theta)),
                  (stan::math::poisson_binomial_log<false>(y, theta)));
}

TEST(ProbDistributionsPoissonBinomial, log_matches_lpmf_vectorial) {
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  int y = 2;
  vec theta(3, 1);
  theta << 0.5, 0.2, 0.7;
  std::vector<vec> thetas{theta, theta};

  EXPECT_FLOAT_EQ((stan::math::poisson_binomial_lpmf(y, thetas)),
                  (stan::math::poisson_binomial_log(y, thetas)));

  EXPECT_FLOAT_EQ((stan::math::poisson_binomial_lpmf<true>(y, thetas)),
                  (stan::math::poisson_binomial_log<true>(y, thetas)));

  EXPECT_FLOAT_EQ((stan::math::poisson_binomial_lpmf<false>(y, thetas)),
                  (stan::math::poisson_binomial_log<false>(y, thetas)));
}
