#include <stan/math/rev.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <boost/math/differentiation/finite_difference.hpp>

TEST(mathMixScalFun, neg_binomial_lpmf_derivatives) {
  auto f1 = [](const auto& alpha, const auto& beta) {
    return stan::math::neg_binomial_lpmf(0, alpha, beta);
  };
  auto f2 = [](const auto& alpha, const auto& beta) {
    return stan::math::neg_binomial_lpmf(6, alpha, beta);
  };

  stan::test::expect_ad(f1, 1.5, 4.1);
  stan::test::expect_ad(f1, std::vector<double>({1.2, 2.0}),
                        std::vector<double>({1.1, 1.2}));
  stan::test::expect_ad(f2, 1.5, 4.1);
  stan::test::expect_ad(f2, std::vector<double>({1.7, 2.0}),
                        std::vector<double>({1.1, 2.3}));
}
