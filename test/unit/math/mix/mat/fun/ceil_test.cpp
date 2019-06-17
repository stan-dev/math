#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, ceil) {
  auto f = [](const auto& x1) { return stan::math::ceil(x1); };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2.1, -0.2, 1.1, 1.51, 3.1);
}
