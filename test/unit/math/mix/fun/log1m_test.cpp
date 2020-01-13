#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, log1m) {
  auto f = [](const auto& x1) { return stan::math::log1m(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -21.5, -21, -1, 0.9, 3);
}
