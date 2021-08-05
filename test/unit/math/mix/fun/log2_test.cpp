#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, log2) {
  auto f = [](const auto& x1) { return stan::math::log2(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -0.2, 1e-3, 1, 1.3, 3, 3.7, 10, 10.2,
                                      1e6);
}
