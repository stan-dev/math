#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

TEST(ProbPoissonLog, cdfs_match_poisson_at_exp_alpha) {
  using stan::math::poisson_cdf;
  using stan::math::poisson_lccdf;
  using stan::math::poisson_lcdf;
  using stan::math::poisson_log_cdf;
  using stan::math::poisson_log_lccdf;
  using stan::math::poisson_log_lcdf;
  std::vector<int> n{0, 3, 17, 82};
  std::vector<double> alpha{-2.0, 0.5, std::log(13.0), std::log(42.0)};
  std::vector<double> lambda;
  for (double a : alpha) {
    lambda.push_back(std::exp(a));
  }
  EXPECT_NEAR(poisson_lcdf(n, lambda), poisson_log_lcdf(n, alpha), 1e-13);
  EXPECT_NEAR(poisson_lccdf(n, lambda), poisson_log_lccdf(n, alpha), 1e-13);
  EXPECT_NEAR(poisson_cdf(n, lambda), poisson_log_cdf(n, alpha), 1e-15);
}

TEST(ProbPoissonLog, cdfs_edge_values) {
  using stan::math::poisson_log_cdf;
  using stan::math::poisson_log_lccdf;
  using stan::math::poisson_log_lcdf;
  const double inf = std::numeric_limits<double>::infinity();
  const int n_max = std::numeric_limits<int>::max();
  // rate exp(alpha) of 0: P(N <= n) = 1 for n >= 0
  EXPECT_EQ(0.0, poisson_log_lcdf(3, -inf));
  EXPECT_EQ(-inf, poisson_log_lccdf(3, -inf));
  EXPECT_EQ(1.0, poisson_log_cdf(3, -inf));
  EXPECT_EQ(-inf, poisson_log_lccdf(3, -800.0));
  // infinite rate, also from overflow of exp(alpha)
  EXPECT_EQ(-inf, poisson_log_lcdf(3, inf));
  EXPECT_EQ(-inf, poisson_log_lcdf(3, 800.0));
  EXPECT_EQ(0.0, poisson_log_lccdf(3, inf));
  EXPECT_EQ(0.0, poisson_log_cdf(3, 800.0));
  // a negative count has P(N <= n) = 0
  EXPECT_EQ(-inf, poisson_log_lcdf(-1, 1.0));
  EXPECT_EQ(0.0, poisson_log_lccdf(-1, 1.0));
  EXPECT_EQ(0.0, poisson_log_cdf(-1, 1.0));
  // n = INT_MAX stands for an infinite count
  EXPECT_EQ(0.0, poisson_log_lcdf(n_max, 1.0));
  EXPECT_EQ(-inf, poisson_log_lccdf(n_max, 1.0));
  EXPECT_EQ(1.0, poisson_log_cdf(n_max, 1.0));
  // an element with probability 1 keeps the other elements
  std::vector<int> n{-1, 3};
  EXPECT_DOUBLE_EQ(poisson_log_lccdf(3, 1.0), poisson_log_lccdf(n, 1.0));
  std::vector<int> m{n_max, 3};
  EXPECT_DOUBLE_EQ(poisson_log_lcdf(3, 1.0), poisson_log_lcdf(m, 1.0));
  EXPECT_DOUBLE_EQ(poisson_log_cdf(3, 1.0), poisson_log_cdf(m, 1.0));
}

TEST(ProbPoissonLog, cdfs_throw_for_nan_alpha) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(stan::math::poisson_log_cdf(3, nan), std::domain_error);
  EXPECT_THROW(stan::math::poisson_log_lcdf(3, nan), std::domain_error);
  EXPECT_THROW(stan::math::poisson_log_lccdf(3, nan), std::domain_error);
}
