#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/ad_tolerances.hpp>

TEST(MathMixMatFun, addDiag) {
  using stan::test::relative_tolerance;
  auto f
      = [](const auto& x, const auto& y) { return stan::math::add_diag(x, y); };

  double d = 1e-10;

  Eigen::MatrixXd m11(1, 1);
  m11 << 42;
  Eigen::VectorXd v1(1);
  v1 << 12;

  Eigen::MatrixXd m00(0, 0);
  Eigen::VectorXd v0(0);

  Eigen::MatrixXd m22(2, 2);
  m22 << 1, 2, 12, 19;
  Eigen::VectorXd v2(2);
  v2 << 0, 1;

  Eigen::MatrixXd m23(2, 3);
  m23 << 1, 1, 1, 1, 1, 1;

  stan::test::ad_tolerances tol;
  tol.hessian_hessian_ = relative_tolerance(1e-4, 1e-3);
  tol.hessian_fvar_hessian_ = relative_tolerance(1e-4, 1e-3);
  // these are OK calls
  stan::test::expect_ad(tol, f, m00, d);
  stan::test::expect_ad(tol, f, m00, v0);
  stan::test::expect_ad(tol, f, m11, d);
  stan::test::expect_ad(tol, f, m11, v1);
  stan::test::expect_ad(tol, f, m22, d);
  stan::test::expect_ad(tol, f, m22, v2);
  stan::test::expect_ad(tol, f, m23, d);
  stan::test::expect_ad(tol, f, m23, v2);

  // these throw
  stan::test::expect_ad(tol, f, m11, v2);
  stan::test::expect_ad(tol, f, m22, v1);
  stan::test::expect_ad(tol, f, m23, v1);
}
