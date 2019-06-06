#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, lgamma) {
  auto f = [](const auto& x1) { return stan::math::lgamma(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
