#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, matrixExpPade) {
  using stan::test::relative_tolerance;
  auto f = [](const auto& x) { return stan::math::matrix_exp_pade(x); };

  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, m00);

  Eigen::MatrixXd m02(0, 2);
  stan::test::expect_ad(f, m02);

  Eigen::MatrixXd a1(1, 1);
  a1 << 1;
  stan::test::expect_ad(f, a1);

  double a = -1;
  double b = -17;
  Eigen::MatrixXd d(2, 2);
  d << -2 * a + 3 * b, 1.5 * a - 1.5 * b, -4 * a + 4 * b, 3 * a - 2 * b;
  stan::test::expect_ad(f, d);

  stan::test::ad_tolerances tols2;
  tols2.hessian_hessian_ = relative_tolerance(5e-4, 1e-3);
  tols2.hessian_fvar_hessian_ = relative_tolerance(5e-4, 1e-3);

  a = -1;
  b = 2;
  double c = 1;
  Eigen::MatrixXd e(3, 3);
  e << -24 * a + 40 * b - 15 * c, 18 * a - 30 * b + 12 * c,
      5 * a - 8 * b + 3 * c, 20 * b - 20 * c, -15 * b + 16 * c, -4 * b + 4 * c,
      -120 * a + 120 * b, 90 * a - 90 * b, 25 * a - 24 * b;
  stan::test::expect_ad(tols2, f, e);

  // replace original random tests
  for (const auto& x : stan::test::ar_test_cov_matrices(1, 3, 0.0)) {
    stan::test::expect_ad(tols2, f, x);
  }
  for (const auto& x : stan::test::ar_test_cov_matrices(1, 3, 0.9)) {
    stan::test::expect_ad(tols2, f, x);
  }
}
