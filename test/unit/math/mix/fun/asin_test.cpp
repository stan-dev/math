#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, asin) {
  auto f = [](const auto& x1) {
    using stan::math::asin;
    return asin(x1);
  };
  // can't autodiff asin through integers
  for (auto x : stan::test::internal::common_nonzero_args())
    stan::test::expect_unary_vectorized(f, x);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -0.2, 0.7, 1.3, 3.4, 5);
  stan::test::expect_complex_common(f);
}
