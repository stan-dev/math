#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, trigamma) {
  auto f = [](const auto& x1) { return stan::math::trigamma(x1); };
  stan::test::expect_unary_vectorized(f, -0.5, 0.5);
  // -0.9, 0, 1, 0.5, 1.3, 5, 10, 19.2);
}
