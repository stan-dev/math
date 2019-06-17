#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, cosh) {
  auto f = [](const auto& x1) { return stan::math::cosh(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -1.2, -0.2, 0.5, 1, 1.3,
                                      1.5);
}
