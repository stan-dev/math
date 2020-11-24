#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, traceQuadForm) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::trace_quad_form(x, y);
  };

  Eigen::MatrixXd a11(1, 1);
  a11 << 1;
  Eigen::VectorXd a1(1);
  a1 << 1;
  stan::test::expect_ad(f, a11, a11);
  stan::test::expect_ad(f, a11, a1);

  stan::test::expect_ad_matvar(f, a11, a11);
  stan::test::expect_ad_matvar(f, a11, a1);

  Eigen::MatrixXd m00(0, 0);
  Eigen::VectorXd v0(0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad(f, m00, v0);

  stan::test::expect_ad_matvar(f, m00, m00);
  stan::test::expect_ad_matvar(f, m00, v0);

  stan::test::ad_tolerances tols;
  tols.hessian_fvar_hessian_ = 1e0;
  tols.hessian_hessian_ = 1e0;

  Eigen::MatrixXd a22(2, 2);
  a22 << 2, 3, 3, 7;
  Eigen::MatrixXd b22(2, 2);
  b22 << 2, 3, 5, 7;
  Eigen::VectorXd a2(2);
  a2 << 2, 3;
  stan::test::expect_ad(tols, f, a22, a22);
  stan::test::expect_ad(tols, f, a22, b22);
  stan::test::expect_ad(tols, f, a22, a2);

  stan::test::expect_ad_matvar(tols, f, a22, a22);
  stan::test::expect_ad_matvar(tols, f, a22, b22);
  stan::test::expect_ad_matvar(tols, f, a22, a2);

  Eigen::MatrixXd a44(4, 4);
  a44 << 2, 3, 4, 5, 6, 10, 2, 2, 7, 2, 7, 1, 8, 2, 1, 112;
  Eigen::MatrixXd b42(4, 2);
  b42 << 100, 10, 0, 1, -3, -3, 5, 2;
  stan::test::expect_ad(tols, f, a44, b42);

  stan::test::expect_ad_matvar(tols, f, a44, b42);

  // size mismatch exceptions
  stan::test::expect_ad(f, a44, b22);
  stan::test::expect_ad(f, a44, a2);

  stan::test::expect_ad_matvar(f, a44, b22);
  stan::test::expect_ad_matvar(f, a44, a2);
}
