#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, expm1) {
  auto f = [](const auto& x1) { return stan::math::expm1(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -0.5, -0.2, 0, 1.0, 1, 1.3,
                                      3);
}
