#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, floor) {
  auto f = [](const auto& x1) { return stan::math::floor(x1); };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2.1 - 0.5, -0.2, 1.1, 1.5,
                                      179.2);
}
