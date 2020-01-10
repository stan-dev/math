#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, trunc) {
  auto f = [](const auto& x1) { return stan::math::trunc(x1); };
  // can't diff trunc through integers
  for (auto x : stan::test::internal::common_nonzero_args())
    stan::test::expect_unary_vectorized(f, x);
  stan::test::expect_unary_vectorized(f, -15.2, 0.5, 1.3, 2.4);
}
