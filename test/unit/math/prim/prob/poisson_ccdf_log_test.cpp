#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbPoisson, ccdf_log_matches_lccdf) {
  int y = 3;
  double lambda = 2.3;

  EXPECT_FLOAT_EQ((stan::math::poisson_lccdf(y, lambda)),
                  (stan::math::poisson_ccdf_log(y, lambda)));
  EXPECT_FLOAT_EQ((stan::math::poisson_lccdf<int, double>(y, lambda)),
                  (stan::math::poisson_ccdf_log<int, double>(y, lambda)));
}

// A negative count has P(N > n) = 1 and contributes 0; the other elements
// of the vector are kept.
TEST(ProbPoisson, lccdf_negative_count_contributes_zero) {
  using stan::math::poisson_lccdf;
  EXPECT_EQ(0.0, poisson_lccdf(-1, 5.0));
  EXPECT_EQ(0.0, poisson_lccdf(-1, 0.0));
  std::vector<int> n{-1, 1};
  EXPECT_DOUBLE_EQ(poisson_lccdf(1, 5.0), poisson_lccdf(n, 5.0));
  std::vector<double> lambda{0.0, 5.0};
  EXPECT_DOUBLE_EQ(poisson_lccdf(1, 5.0), poisson_lccdf(n, lambda));
}
