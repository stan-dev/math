#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::math::poisson_binomial_lccdf;
using vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

TEST(ProbDistributionsPoissonBinomial, lccdf_works_on_scalar_arguments) {
  vec p(3, 1);
  p << 0.5, 0.2, 0.7;

  EXPECT_NEAR(-0.127833, poisson_binomial_lccdf(0, p), 0.001);
  EXPECT_NEAR(-0.798508, poisson_binomial_lccdf(1, p), 0.001);
  EXPECT_NEAR(-2.65926, poisson_binomial_lccdf(2, p), 0.001);
  EXPECT_NEAR(-36.7368, poisson_binomial_lccdf(3, p), 0.001);
}

TEST(ProbDistributionsPoissonBinomial, lccdf_works_on_vectorial_y) {
  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<int> y{2, 2};

  EXPECT_NEAR(-2.65926 * 2, poisson_binomial_lccdf(y, p), 0.001);
}

TEST(ProbDistributionsPoissonBinomial, lccdf_works_on_vectorial_y_and_theta) {
  vec p(3, 1);
  p << 0.5, 0.2, 0.7;
  std::vector<int> y{2, 1};
  std::vector<vec> ps{p, p};

  EXPECT_NEAR(-2.65926 - 0.798508, poisson_binomial_lccdf(y, ps), 0.001);
}
