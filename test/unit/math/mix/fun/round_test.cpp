#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>

TEST_F(mathMix, round) {
  auto f = [](const auto& x1) { return stan::math::round(x1); };
  // can't autodiff round through integers
  for (auto x : stan::test::internal::common_nonzero_args())
    stan::test::expect_unary_vectorized(f, x);
  stan::test::expect_unary_vectorized(f, -2.6, -0.2, 0.4, 1.49, 1.51, 2.4);
}
