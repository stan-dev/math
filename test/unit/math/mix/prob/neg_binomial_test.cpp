#include <stan/math/rev.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <boost/math/differentiation/finite_difference.hpp>

TEST(mathMixScalFun, neg_binomial_lpmf_derivatives) {
  auto f = [](const int y) {
    return [=](const auto& alpha, const auto& beta) {
      return stan::math::neg_binomial_lpmf(y, alpha, beta);
    };
  };

  stan::test::expect_ad(f(0), 1.5, 4.1);
  stan::test::expect_ad(f(0), std::vector<double>({1.2, 2.0}),
                        std::vector<double>({1.1, 1.2}));
  stan::test::expect_ad(f(6), 1.5, 4.1);
  stan::test::expect_ad(f(6), std::vector<double>({1.7, 2.0}),
                        std::vector<double>({1.1, 2.3}));
  stan::test::expect_ad(f(13), 1e11, 1e10);
  stan::test::expect_ad(f(13), std::vector<double>({1e11, 1e10}),
                        std::vector<double>({1e10, 1e9}));
}
