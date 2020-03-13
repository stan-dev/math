#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, quadFormDiag) {
  using stan::test::relative_tolerance;

  auto f = [](const auto& x, const auto& y) {
    return stan::math::quad_form_diag(x, y);
  };
  Eigen::MatrixXd m00(0, 0);
  Eigen::VectorXd v0(0);

  Eigen::MatrixXd m11(1, 1);
  m11 << 1;
  Eigen::VectorXd v1(1);
  v1 << 2;

  Eigen::MatrixXd m22(2, 2);
  m22 << 2, 3, 4, 5;
  Eigen::VectorXd v2(2);
  v2 << 100, 10;

  Eigen::MatrixXd m33(3, 3);
  m33 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Eigen::VectorXd v3(3);
  v3 << 1, 2, 3;

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = relative_tolerance(5e-4, 1e-3);
  tols.hessian_fvar_hessian_ = relative_tolerance(5e-4, 1e-3);

  // matched sizes
  stan::test::expect_ad(tols, f, m00, v0);
  stan::test::expect_ad(tols, f, m11, v1);
  stan::test::expect_ad(tols, f, m22, v2);
  stan::test::expect_ad(tols, f, m33, v3);

  // exceptions from mismatched sizes
  stan::test::expect_ad(tols, f, m33, v2);
  stan::test::expect_ad(tols, f, m22, v3);
}
