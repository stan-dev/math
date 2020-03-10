#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, tan) {
  auto f = [](const auto& x) {
    using stan::math::tan;
    return tan(x);
  };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2, -0.5, 0.5, 1.5, 3, 4.4);
  stan::test::expect_complex_common(f);
}
