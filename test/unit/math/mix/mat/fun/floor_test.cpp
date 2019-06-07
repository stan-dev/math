#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, floor) {
  auto f = [](const auto& x1) { return stan::math::floor(x1); };
  stan::test::expect_common_nonzero_unary_vectorized(f);
}
