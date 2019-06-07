#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, tan) {
  auto f = [](const auto& x1) { return stan::math::tan(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
