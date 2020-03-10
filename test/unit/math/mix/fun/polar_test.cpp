#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <iostream>

TEST(mathMixMatFun, polar) {
  auto f = [](const auto& r, const auto& theta) {
    using stan::math::polar;
    return polar(r, theta);
  };
  stan::test::expect_ad(f, 0.5, 0.5);
  // TODO(carpente): these need to work
  // stan::test::expect_ad(f, 1, 0.5);
  // stan::test::expect_ad(f, 0.5, 1);
  // stan::test::expect_ad(f, 1, 1);
  // stan::test::expect_common_binary(f);

  // need to fix this, too, because it's failing
  // EXPECT_TRUE(stan::is_stan_scalar<std::complex<double>>::value);
}
