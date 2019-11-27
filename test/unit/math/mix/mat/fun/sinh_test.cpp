#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, sinh) {
  auto f = [](const auto& x1) { return stan::math::sinh(x1); };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2, -1.2, -0.5, -0.2, 0.5, 1.3, 1.5,
                                      3);
}
