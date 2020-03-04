#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <vector>

TEST(mixScalFun, abs) {
  auto f = [](const auto& x) { return abs(x); };
  stan::test::expect_common_nonzero_unary(f);
  stan::test::expect_value(f, 0);
  stan::test::expect_value(f, 0.0);

  stan::test::expect_ad(f, -3);
  stan::test::expect_ad(f, -2);
  stan::test::expect_ad(f, 2);

  stan::test::expect_ad(f, -17.3);
  stan::test::expect_ad(f, -0.68);
  stan::test::expect_ad(f, 0.68);
  stan::test::expect_ad(f, 2.0);
  stan::test::expect_ad(f, 4.0);

  // not differentiable at zero
  for (double re :
       std::vector<double>{-4, -2.5, -1.5, -0.3, -0.0, 0.0, 1.3, 2.1, 3.9}) {
    for (double im :
         std::vector<double>{-4, -2.5, -1.5, -0.3, -0.0, 0.0, 1.3, 2.1, 3.9}) {
      if ((re == 0.0 || re == -0.0) && (im == 0.0 || im == -0.0)) {
        continue;
      }
      stan::test::expect_ad(f, std::complex<double>(re, im));
    }
  }
}
