#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, Phi) {
  auto f = [](const auto& x1) { return stan::math::Phi(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -27.5, 27.5);
  for (double x = -37.5; x <= 10; x += 0.5)
    stan::test::expect_unary_vectorized(x);
}
