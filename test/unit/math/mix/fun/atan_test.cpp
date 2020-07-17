#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <vector>

TEST(mathMixMatFun, atan) {
  auto f = [](const auto& x) {
    using stan::math::atan;
    return atan(x);
  };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -2, -0.2, 0.5, 1, 1.3, 1.5, 3);
  // avoid 0 imaginary component where autodiff doesn't work
  for (double re : std::vector<double>{-0.2, 0, 0.3}) {
    for (double im : std::vector<double>{-0.3, 0.2}) {
      stan::test::expect_ad(f, std::complex<double>{re, im});
    }
  }
}
