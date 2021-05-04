#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/ad_tolerances.hpp>

TEST(MathMixMatFun, quadForm) {
  using stan::test::relative_tolerance;
  auto f = [](const auto& x, const auto& y) {
    return stan::math::quad_form(x, y);
  };

  Eigen::MatrixXd a00;
  Eigen::MatrixXd a02(0, 2);
  Eigen::MatrixXd a11(1, 1);
  a11 << 1;
  Eigen::MatrixXd b11(1, 1);
  b11 << -2;
  Eigen::MatrixXd a22(2, 2);
  a22 << 1, 2, 3, 4;
  Eigen::MatrixXd b22(2, 2);
  b22 << -3, -2, -10, 112;
  Eigen::MatrixXd b23(2, 3);
  b23 << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd b42(4, 2);
  b42 << 100, 10, 0, 1, -3, -3, 5, 2;
  Eigen::MatrixXd a44(4, 4);
  a44 << 2, 3, 4, 5, 6, 10, 2, 2, 7, 2, 7, 1, 8, 2, 1, 112;

  Eigen::VectorXd v0(0);
  Eigen::VectorXd v1(1);
  v1 << 42;
  Eigen::VectorXd v2(2);
  v2 << -3, 13;
  Eigen::VectorXd v4(4);
  v4 << 100, 0, -3, 5;

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 2e-1;
  tols.hessian_fvar_hessian_ = 2e-1;

  stan::test::expect_ad(f, a00, a00);
  stan::test::expect_ad(f, a00, a02);
  stan::test::expect_ad(f, a02, a22);

  stan::test::expect_ad_matvar(f, a00, a00);
  stan::test::expect_ad_matvar(f, a00, a02);
  stan::test::expect_ad_matvar(f, a02, a22);

  stan::test::expect_ad(f, a11, b11);
  stan::test::expect_ad(tols, f, a22, b22);
  stan::test::expect_ad(f, a22, b23);
  stan::test::expect_ad(tols, f, a44, b42);

  stan::test::expect_ad_matvar(f, a11, b11);
  stan::test::expect_ad_matvar(tols, f, a22, b22);
  stan::test::expect_ad_matvar(f, a22, b23);
  stan::test::expect_ad_matvar(tols, f, a44, b42);

  stan::test::expect_ad(f, a00, v0);
  stan::test::expect_ad(f, a11, v1);
  stan::test::expect_ad(f, a22, v2);
  stan::test::expect_ad(tols, f, a44, v4);

  stan::test::expect_ad_matvar(f, a00, v0);
  stan::test::expect_ad_matvar(f, a11, v1);
  stan::test::expect_ad_matvar(f, a22, v2);
  stan::test::expect_ad_matvar(tols, f, a44, v4);
}
