#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, acosh) {
  auto f = [](const auto& x1) { return stan::math::acosh(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, 1.5, 1, 3.2, 5, 10, 12.9);
}
