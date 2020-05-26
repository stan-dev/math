#include <test/unit/math/test_ad.hpp>
#include <stan/math/mix.hpp>
#include <gtest/gtest.h>

TEST(mathMixMatFun, lambert_w0) {
  auto f = [](const auto& x1) {
    using stan::math::lambert_w0;
    return lambert_w0(x1);
  };
  stan::test::expect_ad_vectorized(f, -0.3);
  stan::test::expect_ad_vectorized(f, -0.1);
  stan::test::expect_ad_vectorized(f, 0.0);
  stan::test::expect_ad_vectorized(f, 1);
  stan::test::expect_ad_vectorized(f, 10);
  stan::test::expect_ad_vectorized(f, 20);
  stan::math::recover_memory();
}
