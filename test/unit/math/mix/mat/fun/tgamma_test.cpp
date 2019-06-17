#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, tgamma) {
  auto f = [](const auto& x1) { return stan::math::tgamma(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -1.2, -0.5, 0.5, 1.5, 3.5);
}
