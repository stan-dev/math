#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, floor) {
  auto f = [](const auto& x1) { return stan::math::floor(x1); };
  // can't autodiff floor through integers
  for (auto x : stan::test::internal::common_nonzero_args())
    stan::test::expect_unary_vectorized(f, x);
  stan::test::expect_unary_vectorized(f, -2.6, -2.1 - 0.5, -0.2, 1.1, 1.5,
                                      179.2);
}
