#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <vector>

namespace {

// largest |r_j| over the integers j in [lo, hi], by evaluation of every
// ratio r_j = z (a1 + j)(a2 + j)(a3 + j) / ((b1 + j)(b2 + j)(1 + j))
double max_abs_ratio(const std::array<double, 3>& a,
                     const std::array<double, 2>& b, double z, int lo, int hi) {
  double max_r = 0.0;
  for (int j = lo; j <= hi; ++j) {
    const double r = z * (a[0] + j) * (a[1] + j) * (a[2] + j)
                     / ((b[0] + j) * (b[1] + j) * (1.0 + j));
    max_r = std::fmax(max_r, std::fabs(r));
  }
  return max_r;
}

}  // namespace

TEST(MathPrimScalFun, hypergeometric_3F2_min_abs) {
  using stan::math::internal::hypergeometric_3F2_min_abs;
  // no zero of p + j in [1, 10]
  EXPECT_DOUBLE_EQ(3.5, hypergeometric_3F2_min_abs(2.5, 1, 10));
  EXPECT_DOUBLE_EQ(5.25, hypergeometric_3F2_min_abs(-15.25, 1, 10));
  // zero at j = 5.25, between the integers 5 and 6
  EXPECT_DOUBLE_EQ(0.25, hypergeometric_3F2_min_abs(-5.25, 1, 10));
  EXPECT_DOUBLE_EQ(0.0, hypergeometric_3F2_min_abs(-3.0, 1, 10));
}

TEST(MathPrimScalFun, hypergeometric_3F2_ratio_bound_is_upper_bound) {
  using stan::math::internal::hypergeometric_3F2_ratio_bound;
  struct series {
    std::array<double, 3> a;
    std::array<double, 2> b;
    double z;
    int lo;
    int hi;
  };
  const std::vector<series> cases = {
      // numerator and denominator parameters with zeros inside [lo, hi]
      {{0.5, -2.5, -10.0}, {-3.5, 1.5}, 0.7, 0, 9},
      {{2.0, 3.0, -20.0}, {-7.25, 0.5}, -0.9, 0, 19},
      {{-4.5, 7.0, -30.0}, {-12.5, -40.5}, 1.0, 3, 29},
      // a tiny first term before large terms
      {{1e-23, 300.0, -100.0}, {1.0, -150.5}, 1.0, 0, 99},
      {{1e-23, 300.0, -100.0}, {1.0, -150.5}, 1.0, 60, 99},
      // the series of beta_binomial_lcdf(3 | 10, 0.5, 0.5)
      {{1.0, 4.5, -6.0}, {5.0, -5.5}, 1.0, 0, 5},
  };
  for (const auto& s : cases) {
    const double bound
        = hypergeometric_3F2_ratio_bound(s.a, s.b, std::fabs(s.z), s.lo, s.hi);
    EXPECT_GE(bound * (1 + 1e-14), max_abs_ratio(s.a, s.b, s.z, s.lo, s.hi));
  }
}

TEST(MathPrimScalFun, hypergeometric_3F2_ratio_bound_is_tight) {
  using stan::math::internal::hypergeometric_3F2_ratio_bound;
  // the series of beta_binomial_lcdf(5000 | 10000, 1000, 3000): every factor
  // decreases in j, so the bound equals the ratio at j = lo
  const std::array<double, 3> a{1.0, 6001.0, -4999.0};
  const std::array<double, 2> b{5002.0, -7998.0};
  for (int lo : {0, 100, 1000, 4000}) {
    const double max_r = max_abs_ratio(a, b, 1.0, lo, 4998);
    EXPECT_NEAR(max_r, hypergeometric_3F2_ratio_bound(a, b, 1.0, lo, 4998),
                1e-14 * max_r);
  }
}

TEST(MathPrimScalFun, hypergeometric_3F2_tail_sums) {
  using stan::math::internal::hypergeometric_3F2_tail_sums;
  for (double rho : {0.0, 0.3, 0.9, 1.0, 1.7}) {
    for (int n : {1, 5, 40}) {
      double sum1 = 0.0;
      double sum2 = 0.0;
      for (int d = 1; d <= n; ++d) {
        sum1 += std::pow(rho, d);
        sum2 += d * std::pow(rho, d);
      }
      const auto sums = hypergeometric_3F2_tail_sums(rho, n);
      EXPECT_GE(sums.first, sum1 * (1 - 1e-15));
      EXPECT_GE(sums.second, sum2 * (1 - 1e-15));
    }
  }
  EXPECT_EQ(0.0, hypergeometric_3F2_tail_sums(0.5, 0).first);
  EXPECT_EQ(0.0, hypergeometric_3F2_tail_sums(0.5, 0).second);
}
