#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, exp) {
  auto f = [](const auto& x1) { return stan::math::exp(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -15.2, -10, -0.5, 0.5, 1, 1.0, 1.3, 5,
                                      10);
}
