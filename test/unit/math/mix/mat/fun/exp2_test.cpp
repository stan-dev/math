#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, exp2) {
  auto f = [](const auto& x1) { return stan::math::exp2(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
