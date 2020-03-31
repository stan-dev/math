#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixMatFun, sqrt) {
  auto f = [](const auto& x1) {
    using stan::math::sqrt;
    return sqrt(x1);
  };
  stan::test::expect_common_nonzero_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -6, -5.2, 1.3, 7, 10.7, 36, 1e6);

  // undefined with 0 in denominator
  stan::test::expect_ad(f, std::complex<double>(0.9, 0.8));
  for (double im : std::vector<double>{-1.3, 2.3}) {
    for (double re : std::vector<double>{-3.6, -0.0, 0.0, 0.5}) {
      stan::test::expect_ad(f, std::complex<double>(re, im));
    }
  }
}
