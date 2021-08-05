#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, logit) {
  auto f = [](const auto& x1) { return stan::math::logit(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -1.2, -0.5, 0.01, 0.5, 0.99, 1.5);
}
