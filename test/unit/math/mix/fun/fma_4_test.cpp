#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, fma_matrix) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    return stan::math::fma(x1, x2, x3).eval();
  };

  double xd = 1.0;
  Eigen::MatrixXd xm(2, 2);
  xm << 1.0, 2.0, -1.0, 1.1;

  double yd = 2.0;
  Eigen::MatrixXd ym(2, 2);
  ym << 1.0, 2.0, -1.0, 1.1;

  double zd = 3.0;
  Eigen::MatrixXd zm(2, 2);
  zm << 1.0, 2.0, -1.0, 1.1;

  stan::test::expect_ad(f, xd, yd, zm);
  stan::test::expect_ad(f, xd, ym, zd);
  stan::test::expect_ad(f, xd, ym, zm);
  stan::test::expect_ad(f, xm, yd, zd);
  stan::test::expect_ad(f, xm, yd, zm);
  stan::test::expect_ad(f, xm, ym, zd);

  stan::test::expect_ad(f, xm, ym, zm);

  stan::test::expect_ad_matvar(f, xd, yd, zm);
  stan::test::expect_ad_matvar(f, xd, ym, zd);
  stan::test::expect_ad_matvar(f, xd, ym, zm);
  stan::test::expect_ad_matvar(f, xm, yd, zd);
  stan::test::expect_ad_matvar(f, xm, yd, zm);
  stan::test::expect_ad_matvar(f, xm, ym, zd);
  stan::test::expect_ad_matvar(f, xm, ym, zm);
}
