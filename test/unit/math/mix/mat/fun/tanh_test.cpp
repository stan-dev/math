#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, tanh) {
  auto f = [](const auto& x1) { return stan::math::tanh(x1); };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -1.2, -0.5, 0.5, 1.5);
}
