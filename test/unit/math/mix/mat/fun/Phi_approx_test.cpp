#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, PhiApprox) {
  auto f = [](const auto& x1) { return stan::math::Phi_approx(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -3.0, 1, 1.3, 3);
}
