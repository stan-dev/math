#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>


TEST(ProbDistributionsPoissonBinomial, lcdf_works_on_scalar_arguments) {
  using stan::math::poisson_binomial_lcdf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;

  EXPECT_NEAR(-2.12026, poisson_binomial_lcdf(0, p), 0.001);
  EXPECT_NEAR(-0.597837, poisson_binomial_lcdf(1, p), 0.001);
  EXPECT_NEAR(-0.0725707, poisson_binomial_lcdf(2, p), 0.001);
  EXPECT_NEAR(0, poisson_binomial_lcdf(3, p), 0.001);
}

TEST(ProbDistributionsPoissonBinomial, lcdf_works_on_vectorial_y) {
  using stan::math::poisson_binomial_lcdf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<int> y{2, 2};

  EXPECT_NEAR(-0.0725707 * 2, poisson_binomial_lcdf(y, p), 0.001);
}

TEST(ProbDistributionsPoissonBinomial, lcdf_works_on_vectorial_y_and_theta) {
  using stan::math::poisson_binomial_lcdf;
  using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;

  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<int> y{2, 1};
  std::vector<vec> ps{p, p};

  EXPECT_NEAR(-0.0725707 -0.597837, poisson_binomial_lcdf(y, ps), 0.001);
}
