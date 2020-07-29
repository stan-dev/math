#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>

TEST(mathMixScalFun, frechet_lpmf_derivatives) {
  auto f = [](const auto& y, const auto& alpha, const auto& beta) {
    return stan::math::frechet_lpdf<false>(y, alpha, beta);
  };

  stan::test::expect_ad(f, 2.0, 1.0, 1.0);
}
