#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, round) {
  auto f = [](const auto& x1) { return stan::math::round(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
