#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, asinh) {
  auto f = [](const auto& x1) { return stan::math::asinh(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -1.2, -0.2, 0.5, 2, -1.2);
}
