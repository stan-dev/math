#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixMatFun, acos) {
  using stan::test::expect_unary_vectorized;
  auto f = [](const auto& x1) { return stan::math::acos(x1); };
  // can't autodiff acos through integers
  for (auto x : stan::test::internal::common_nonzero_args())
    stan::test::expect_unary_vectorized(f, x);
  expect_unary_vectorized(f, -2.2, -0.8, 0.5,
                          1 + std::numeric_limits<double>::epsilon(), 1.5, 3,
                          3.4, 4);
}
