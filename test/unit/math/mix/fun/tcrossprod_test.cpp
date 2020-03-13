#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, tcrossprod) {
  using stan::test::relative_tolerance;
  auto f = [](const auto& y) { return stan::math::tcrossprod(y); };

  Eigen::MatrixXd a00(0, 0);
  stan::test::expect_ad(f, a00);

  Eigen::MatrixXd a11(1, 1);
  a11 << 3;
  stan::test::expect_ad(f, a11);

  Eigen::MatrixXd a22(2, 2);
  a22 << 3, 0, 4, -3;
  stan::test::expect_ad(f, a22);

  Eigen::MatrixXd b22(2, 2);
  b22 << 1, 0, 2, 3;
  stan::test::expect_ad(f, b22);

  Eigen::MatrixXd a13(1, 3);
  a13 << 1, 2, 3;
  stan::test::expect_ad(f, a13);

  Eigen::MatrixXd a31 = a13.transpose();
  stan::test::expect_ad(f, a31);

  Eigen::MatrixXd a23(2, 3);
  a23 << 1, 2, 3, -1, 4, -9;
  stan::test::expect_ad(f, a23);

  Eigen::MatrixXd a32 = a23.transpose();
  stan::test::expect_ad(f, a32);

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = relative_tolerance(2e-3, 2e-4);
  tols.hessian_fvar_hessian_ = relative_tolerance(2e-3, 2e-4);

  Eigen::MatrixXd a33(3, 3);
  a33 << 1, 0, 0, 2, 3, 0, 4, 5, 6;
  stan::test::expect_ad(tols, f, a33);

  Eigen::MatrixXd b33(3, 3);
  b33 << 1, 2, 3, 1, 4, 9, 1, 8, 27;
  stan::test::expect_ad(tols, f, b33);
}
