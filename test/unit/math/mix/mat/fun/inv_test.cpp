#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, inv) {
  auto f = [](const auto& x1) { return stan::math::inv(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
