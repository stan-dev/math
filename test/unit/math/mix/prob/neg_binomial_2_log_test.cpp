#include <stan/math/rev.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <boost/math/differentiation/finite_difference.hpp>

TEST(mathMixScalFun, neg_binomial_2_log_lpmf_derivatives) {
  auto f1 = [](const auto& eta, const auto& phi) {
    return stan::math::neg_binomial_2_log_lpmf(0, eta, phi);
  };
  auto f2 = [](const auto& eta, const auto& phi) {
    return stan::math::neg_binomial_2_log_lpmf(6, eta, phi);
  };

  stan::test::expect_ad(f1, -1.5, 4.1);
  stan::test::expect_ad(f1, 2.0, 1.1);
  stan::test::expect_ad(f2, -1.5, 4.1);
  stan::test::expect_ad(f2, 2.0, 1.1);
}
