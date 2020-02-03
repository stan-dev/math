#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, atanh) {
  auto f = [](const auto& x1) { return stan::math::atanh(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -0.9, 0.5);
}
