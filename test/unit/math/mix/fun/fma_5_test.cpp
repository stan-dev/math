#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, fma_matrix_error) {
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    return stan::math::fma(x1, x2, x3);
  };

  double xd = 1.0;
  Eigen::RowVectorXd xr(2);
  xr << 1.0, 2.0;

  double yd = 2.0;
  Eigen::VectorXd yv(2);
  yv << 2.0, -3.0;
  Eigen::RowVectorXd yr(2);
  yr << 2.0, -3.0;

  double zd = 3.0;
  Eigen::VectorXd zv(2);
  zv << -3.0, 4.0;

  stan::test::expect_ad(f, xd, yr, zv);
  stan::test::expect_ad(f, xr, yd, zv);
  stan::test::expect_ad(f, xr, yv, zd);

  stan::test::expect_ad_matvar(f, xd, yr, zv);
  stan::test::expect_ad_matvar(f, xr, yd, zv);
  stan::test::expect_ad_matvar(f, xr, yv, zd);
}
