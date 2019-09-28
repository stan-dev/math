#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, acosh) {
  auto f = [](const auto& x1) { return stan::math::acosh(x1); };
  for (double x : stan::test::internal::common_args())
    stan::test::expect_unary_vectorized(x);
  stan::test::expect_unary_vectorized(f, 1.5, 3.2, 5, 10, 12.9);
}
