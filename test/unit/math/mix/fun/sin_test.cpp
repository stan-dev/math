#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, sin) {
  auto f = [](const auto& x1) {
    using stan::math::sin;
    return sin(x1);
  };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -0.2, 3, 5, 5.3);
  stan::test::expect_complex_common(f);
}
