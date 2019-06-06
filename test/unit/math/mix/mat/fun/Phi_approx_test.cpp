#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, PhiApprox) {
  auto f = [](const auto& x1) { return stan::math::Phi_approx(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
