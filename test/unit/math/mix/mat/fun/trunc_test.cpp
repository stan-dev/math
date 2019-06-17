#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, trunc) {
  auto f = [](const auto& x1) { return stan::math::trunc(x1); };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -15.2, 0.5, 1.3, 2.4);
}
