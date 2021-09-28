#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, fma_row_vector) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    return stan::math::fma(x1, x2, x3);
  };

  double xd = 1.0;
  Eigen::RowVectorXd xr(2);
  xr << 1.0, 2.0;

  double yd = 2.0;
  Eigen::RowVectorXd yr(2);
  yr << 2.0, -3.0;

  double zd = 3.0;
  Eigen::RowVectorXd zr(2);
  zr << -1.0, 2.0;

  stan::test::expect_ad(f, xd, yd, zr);
  stan::test::expect_ad(f, xd, yr, zd);
  stan::test::expect_ad(f, xd, yr, zr);
  stan::test::expect_ad(f, xr, yd, zd);
  stan::test::expect_ad(f, xr, yd, zr);
  stan::test::expect_ad(f, xr, yr, zd);
  stan::test::expect_ad(f, xr, yr, zr);

  stan::test::expect_ad_matvar(f, xd, yd, zr);
  stan::test::expect_ad_matvar(f, xd, yr, zd);
  stan::test::expect_ad_matvar(f, xd, yr, zr);
  stan::test::expect_ad_matvar(f, xr, yd, zd);
  stan::test::expect_ad_matvar(f, xr, yd, zr);
  stan::test::expect_ad_matvar(f, xr, yr, zd);
  stan::test::expect_ad_matvar(f, xr, yr, zr);
}
