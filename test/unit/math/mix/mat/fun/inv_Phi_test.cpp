#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, invPhi) {
  auto f = [](const auto& x1) { return stan::math::inv_Phi(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
