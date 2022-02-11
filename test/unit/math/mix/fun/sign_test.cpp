#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, sign) {
  auto f = [](const auto& x) { return stan::math::sign(x); };

  stan::test::expect_ad(f, -13.2);
  stan::test::expect_ad(f, -0.1);
  stan::test::expect_ad(f, 0.01);
  stan::test::expect_ad(f, 0.1);
  stan::test::expect_ad(f, 1.5);
  stan::test::expect_ad(f, 100.1);
  stan::test::expect_ad(f, std::numeric_limits<double>::quiet_NaN());

  for (auto x : stan::test::internal::common_nonzero_args())
    stan::test::expect_unary_vectorized(f, x);
}
