#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, cosh) {
  auto f = [](const auto& x1) { return stan::math::cosh(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
