#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, exp2) {
  auto f = [](const auto& x1) { return stan::math::exp2(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -15.2, -10, 1, 1.3, 5, 10);
}
