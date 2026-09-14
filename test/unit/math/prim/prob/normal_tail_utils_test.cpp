#include <stan/math/prim/prob/normal_tail_utils.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>

namespace {

void expect_stable_value(double actual, double expected) {
  EXPECT_NEAR(actual, expected, 1e-12 * std::max(1.0, std::abs(expected)));
}

void expect_normal_tail_terms(double x, double expected_log_cdf,
                              double expected_mills_ratio) {
  const auto terms = stan::math::internal::std_normal_lcdf_and_mills(x);
  expect_stable_value(terms.log_cdf, expected_log_cdf);
  expect_stable_value(terms.mills_ratio, expected_mills_ratio);
}

}  // namespace

TEST(ProbNormalTailUtils, central_and_lower_tail) {
  expect_normal_tail_terms(0.0, -0.69314718055994529,
                           0.79788456080286529);
  expect_normal_tail_terms(-7.0710678118654755, -27.894036726097383,
                           7.2073273273945926);
  expect_normal_tail_terms(-28.284271247461902, -404.26249051466425,
                           28.319538745551821);
  expect_normal_tail_terms(-40.0, -804.6084420137538, 40.024968847206338);
}
