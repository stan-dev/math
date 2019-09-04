#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, sqrt) {
  auto f = [](const auto& x1) { return stan::math::sqrt(x1); };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -6, -5.2, 1.3, 7, 10.7, 36, 1e6);
}
