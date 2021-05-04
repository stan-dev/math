#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, scaleMatrixExpMultiply) {
  using stan::test::relative_tolerance;
  auto f = [](const auto& t, const auto& a, const auto& b) {
    return stan::math::scale_matrix_exp_multiply(t, a, b);
  };

  double t = 1.3;

  // 0 x 0
  Eigen::MatrixXd a00(0, 0);
  Eigen::MatrixXd b00(0, 0);
  stan::test::expect_ad(f, t, a00, b00);
  Eigen::MatrixXd b03(0, 3);
  stan::test::expect_ad(f, t, a00, b03);

  // 1 x 1
  Eigen::MatrixXd a11(1, 1);
  a11 << 0.5;
  Eigen::MatrixXd b11(1, 1);
  b11 << 1.3;
  stan::test::expect_ad(f, t, a11, b11);

  // 1 x 4
  Eigen::MatrixXd b14(1, 4);
  b14 << 1, 2, 3, -1;
  stan::test::expect_ad(f, t, a11, b14);

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = relative_tolerance(1e-4, 1e-2);
  tols.hessian_fvar_hessian_ = relative_tolerance(1e-4, 1e-2);

  // 3 x 1
  Eigen::MatrixXd a33(3, 3);
  a33 << 0.1, 0.4, 0.6, 0.3, 0.3, 0.4, 0.5, 0.1, 0.4;
  Eigen::MatrixXd b31(3, 1);
  b31 << 0.3, 0.5, -0.2;
  stan::test::expect_ad(tols, f, t, a33, b31);

  // 3 x 3
  Eigen::MatrixXd b33(3, 3);
  b33 << -5, 4, -3, 2, 1, 0, -1, 2, -3;
  stan::test::expect_ad(tols, f, t, a33, b33);
}
