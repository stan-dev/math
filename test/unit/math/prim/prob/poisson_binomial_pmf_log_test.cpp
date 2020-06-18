#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsPoissonBinomial, lpmf_works_on_scalar_arguments) {
  using stan::math::poisson_binomial_lpmf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;

  EXPECT_NEAR(-2.12026, poisson_binomial_lpmf(0, p), 0.001);
  EXPECT_NEAR(-0.84397, poisson_binomial_lpmf(1, p), 0.001);
  EXPECT_NEAR(-0.967584, poisson_binomial_lpmf(2, p), 0.001);
  EXPECT_NEAR(-2.65926, poisson_binomial_lpmf(3, p), 0.001);
}

TEST(ProbDistributionsPoissonBinomial, lpmf_works_on_vectorial_y) {
  using stan::math::poisson_binomial_lpmf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<int> y{2, 2};

  EXPECT_NEAR(-0.967584 * 2, poisson_binomial_lpmf(y, p), 0.001);
}

TEST(ProbDistributionsPoissonBinomial, lpmf_works_on_vectorial_y_and_theta) {
  using stan::math::poisson_binomial_lpmf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<int> y{2, 0};
  std::vector<vec> ps{p, p};

  EXPECT_NEAR(-0.967584 - 2.12026, poisson_binomial_lpmf(y, ps), 0.001);
}
