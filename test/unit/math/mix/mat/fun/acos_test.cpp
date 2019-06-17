#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixMatFun, acos) {
  using stan::test::expect_unary_vectorized;
  auto f = [](const auto& x1) { return stan::math::acos(x1); };
  stan::test::expect_common_unary_vectorized(f);
  expect_unary_vectorized(f, -2.2, -1, -0.8, 0.5, 1,
                          1 + std::numeric_limits<double>::epsilon(), 1.5, 3,
                          3.4, 4);
}
