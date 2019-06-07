#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, log10) {
  auto f = [](const auto& x1) { return stan::math::log10(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
