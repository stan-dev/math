#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, invPhi) {
  auto f = [](const auto& x1) { return stan::math::inv_Phi(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, 0.02425, 0.97575);  // breakpoints
  stan::test::expect_unary_vectorized(f, -100.25, -2, 0.01, 0.1, 0.98, 0.5,
                                      2.0);
}
