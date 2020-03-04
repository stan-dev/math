#include <test/unit/math/test_ad.hpp>

TEST(mathMixFun, asinh) {
  auto f = [](const auto& x1) { return stan::math::asinh(x1); };
  //  stan::test::expect_common_unary_vectorized(f);
  // stan::test::expect_unary_vectorized(f, -2.6, -1.2, -0.2, 0.5, 2, -1.2);
  // stan::test::expect_complex_common(f);
}

TEST(mathMixFun, asinhComplex) {
  auto f = [](const auto& z) { return asinh(z); };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;
  // for (double re : std::vector<double>{-1.1, 2.9}) {
  //   for (double im : std::vector<double>{-1.1, 2.9}) {
  //     stan::test::expect_ad(tols, f, std::complex<double>{re, im});
  //   }
  // }
  // stan::test::expect_ad(tols, f, std::complex<double>{1.1, -1.3});
  // stan::test::expect_ad(tols, f, std::complex<double>{-1.1, -1.3});
  stan::test::expect_ad(tols, f, std::complex<double>{0.1, 0.1});

  stan::test::expect_complex_common(f);
}
