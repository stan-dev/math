#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, erfc) {
  auto f = [](const auto& x1) { return stan::math::erfc(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
