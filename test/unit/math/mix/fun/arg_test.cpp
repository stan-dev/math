#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <vector>

TEST(mixScalFun, arg) {
  auto f = [](const auto& x) { return arg(x); };

  // undefined with 0 in denominator
  stan::test::expect_ad(f, std::complex<double>(0.9, 0.8));
  for (double re : std::vector<double>{-3.6, -0.0, 0.0, 0.5}) {
    for (double im : std::vector<double>{-1.3, 2.3}) {
      stan::test::expect_ad(f, std::complex<double>(re, im));
    }
  }
}
