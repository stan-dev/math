#include <test/unit/math/test_ad.hpp>
#include <stan/math/prim/core/complex_base.hpp>

#include <vector>

TEST(mixScalFun, arg) {
  auto f = [](const auto& x) { return arg(x); };

  // undefined with 0 in denominator
  stan::test::expect_ad(f, stan::math::complex<double>(0.9, 0.8));
  for (double re : std::vector<double>{-3.6, -0.0, 0.0, 0.5}) {
    for (double im : std::vectostan::math::complex.3, 2.3}) {
      stan::test::expect_ad(f, stan::math::complex<double>(re, im));
    }
  }
}
