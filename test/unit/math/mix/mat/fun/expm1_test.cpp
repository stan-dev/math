#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, expm1) {
  auto f = [](const auto& x1) { return stan::math::expm1(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
