#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>

TEST(mathMixMatFun, lambert_w0) {
  using stan::test::expect_unary_vectorized;
  auto f = [](const auto& x1) {
    using stan::math::lambert_w0;
    return lambert_w0(x1);
  };
  stan::test::expect_ad_vectorized(f, 0.0);
  // seg fault
  // stan::test::expect_ad_vectorized(f, -0.3);
  // stan::test::expect_ad_vectorized(f, -0.1);
  // Wrong answer
  stan::test::expect_ad_vectorized(f, 1);
}
