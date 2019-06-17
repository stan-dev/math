#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, log1p) {
  auto f = [](const auto& x1) { return stan::math::log1p(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2, -.99, -0.5, 0.5, 7.2, 1000);
}
