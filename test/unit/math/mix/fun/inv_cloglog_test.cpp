#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, invCLogLog) {
  auto f = [](const auto& x1) { return stan::math::inv_cloglog(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -1.2, -0.2, 0.5, 1.3);
}
