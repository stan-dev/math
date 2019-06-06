#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, atan) {
  auto f = [](const auto& x1) { return stan::math::atan(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
