#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, log) {
  auto f = [](const auto& x1) { return stan::math::log(x1); };
  stan::test::expect_common_unary_vectorized(f);
}
