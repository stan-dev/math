#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <limits>
#include <vector>

TEST(mathMixFun, complexPow) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::pow;
    return pow(x1, x2);
  };
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 5e-3;
  tols.hessian_fvar_hessian_ = 5e-3;
  // complex, complex
  for (auto re1 : std::vector<double>{-1.8, 3.4}) {
    for (auto im1 : std::vector<double>{-1.8, 3.4}) {
      for (auto re2 : std::vector<double>{-2.7, 1, 2.3}) {
        for (auto im2 : std::vector<double>{-1.5, 1.2}) {
          stan::test::expect_ad(tols, f, std::complex<double>{re1, im1},
                                std::complex<double>{re2, im2});
        }
      }
    }
  }
  // if real, first arg must be positive
  for (auto re1 : std::vector<double>{3.4}) {
    for (auto re2 : std::vector<double>{-2.7, 1, 2.3}) {
      for (auto im2 : std::vector<double>{-1.5, 1.2}) {
        stan::test::expect_ad(tols, f, re1, std::complex<double>{re2, im2});
      }
    }
  }
  for (auto re1 : std::vector<double>{-1.8, 3.4}) {
    for (auto im1 : std::vector<double>{-1.8, 3.4}) {
      for (auto re2 : std::vector<double>{-2.7, 1, 2.3}) {
        stan::test::expect_ad(tols, f, std::complex<double>{re1, im1}, re2);
      }
    }
  }

  Eigen::Matrix<std::complex<double>, -1, 1> din1(2);
  din1.real() << 0.5, 0.1;
  din1.imag() << 1.6, 5.4;
  Eigen::Matrix<std::complex<double>, -1, 1> din2(2);
  din2.real() << 1.2, 2.3;
  din2.imag() << 8.1, 6.1;

  stan::test::expect_ad_vectorized_binary(tols, f, din1, din2);
}

TEST(mathMixFun, powIntAmbiguityTest) {
  using stan::math::pow;  // included to check ambiguities
  using stan::math::var;
  using std::complex;
  int i = 2;
  double d = 2.5;
  var v = 2.5;
  complex<double> cd = 2.5;
  complex<var> cv = 2.5;

  auto a1 = pow(i, i);
  auto a2 = pow(i, d);
  auto a3 = pow(i, v);
  auto a4 = pow(i, cd);
  auto a5 = pow(i, cv);

  auto b1 = pow(d, i);
  auto b2 = pow(d, d);
  auto b3 = pow(d, v);
  auto b4 = pow(d, cd);
  auto b5 = pow(d, cv);

  auto e1 = pow(v, i);
  auto e2 = pow(v, d);
  auto e3 = pow(v, v);
  auto e4 = pow(v, cd);
  auto e5 = pow(v, cv);

  auto c1 = pow(cd, i);
  auto c2 = pow(cd, d);
  auto c3 = pow(cd, v);
  auto c4 = pow(cd, cd);
  auto c5 = pow(cd, cv);

  auto d1 = pow(cv, i);
  auto d2 = pow(cv, d);
  auto d3 = pow(cv, v);
  auto d4 = pow(cv, cd);
  auto d5 = pow(cv, cv);

  auto e = a1 + a2 + a3 + a4 + a5 + b1 + b2 + b3 + b4 + b5 + c1 + c2 + c3 + c4
           + c5 + d1 + d2 + d3 + d4 + d5;
}

TEST(mathMixFun, powIntAmbiguityTestFvar) {
  using stan::math::fvar;
  using stan::math::pow;  // included to check ambiguities
  using std::complex;
  int i = 2;
  double d = 2.5;
  fvar<double> v = 2.5;
  complex<double> cd = 2.5;
  complex<fvar<double>> cv = 2.5;

  auto a1 = pow(i, i);
  auto a2 = pow(i, d);
  auto a3 = pow(i, v);
  auto a4 = pow(i, cd);
  auto a5 = pow(i, cv);

  auto b1 = pow(d, i);
  auto b2 = pow(d, d);
  auto b3 = pow(d, v);
  auto b4 = pow(d, cd);
  auto b5 = pow(d, cv);

  auto c1 = pow(cd, i);
  auto c2 = pow(cd, d);
  auto c3 = pow(cd, v);
  auto c4 = pow(cd, cd);
  auto c5 = pow(cd, cv);

  auto d1 = pow(cv, i);
  auto d2 = pow(cv, d);
  auto d3 = pow(cv, v);
  auto d4 = pow(cv, cd);
  auto d5 = pow(cv, cv);

  auto e = a1 + a2 + a3 + a4 + a5 + b1 + b2 + b3 + b4 + b5 + c1 + c2 + c3 + c4
           + c5 + d1 + d2 + d3 + d4 + d5;
}
