#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, erf) {
  auto f = [](const auto& x1) { return stan::math::erf(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2, -1, -0.2, 2.6);
}
