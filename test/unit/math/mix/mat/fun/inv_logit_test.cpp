#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, invLogit) {
  auto f = [](const auto& x1) { return stan::math::inv_logit(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
