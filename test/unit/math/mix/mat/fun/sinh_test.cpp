#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, sinh) {
  auto f = [](const auto& x1) { return stan::math::sinh(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
