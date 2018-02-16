#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>

TEST(ProbHypergeometric, log_matches_lpmf) {
  int n = 1;
  int N = 4;
  double a = 1.1;
  double b = 3.5;

  EXPECT_FLOAT_EQ((stan::math::hypergeometric_lpmf(n, N, a, b)),
                  (stan::math::hypergeometric_log(n, N, a, b)));
  EXPECT_FLOAT_EQ((stan::math::hypergeometric_lpmf<true>(n, N, a, b)),
                  (stan::math::hypergeometric_log<true>(n, N, a, b)));
  EXPECT_FLOAT_EQ((stan::math::hypergeometric_lpmf<false>(n, N, a, b)),
                  (stan::math::hypergeometric_log<false>(n, N, a, b)));
  EXPECT_FLOAT_EQ(
      (stan::math::hypergeometric_lpmf<true, int, int, double, double>(n, N, a,
                                                                       b)),
      (stan::math::hypergeometric_log<true, int, int, double, double>(n, N, a,
                                                                      b)));
  EXPECT_FLOAT_EQ(
      (stan::math::hypergeometric_lpmf<false, int, int, double, double>(n, N, a,
                                                                        b)),
      (stan::math::hypergeometric_log<false, int, int, double, double>(n, N, a,
                                                                       b)));
  EXPECT_FLOAT_EQ(
      (stan::math::hypergeometric_lpmf<int, int, double, double>(n, N, a, b)),
      (stan::math::hypergeometric_log<int, int, double, double>(n, N, a, b)));
}
