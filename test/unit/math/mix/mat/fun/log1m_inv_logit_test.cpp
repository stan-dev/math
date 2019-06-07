#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, log1mInvLogit) {
  auto f = [](const auto& x1) { return stan::math::log1m_inv_logit(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
