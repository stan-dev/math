#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, cos) {
  auto f = [](const auto& x1) { return stan::math::cos(x1); };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -0.2, -0.5, 0, 1.5, 3, 5,
                                      5.3);
}
