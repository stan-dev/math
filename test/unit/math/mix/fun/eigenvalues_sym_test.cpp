#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, eigenvaluesSym) {
  using stan::test::relative_tolerance;
  auto f = [](const auto& y) {
    // need to maintain symmetry for finite diffs
    auto a = ((y + y.transpose()) * 0.5).eval();
    return stan::math::eigenvalues_sym(a);
  };

  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, m00);

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 5e-2;
  tols.hessian_fvar_hessian_ = 5e-2;
  // tols.grad_hessian_grad_hessian_ = 5e-1;  // challenging

  // 1 x 1
  Eigen::MatrixXd a11(1, 1);
  a11 << 2.3;
  stan::test::expect_ad(tols, f, a11);

  // 3 x 3
  Eigen::MatrixXd a33(3, 3);
  a33 << 1, 2, 3, 2, 5, 7.9, 3, 7.9, 1.08;
  stan::test::expect_ad(tols, f, a33);

  // asymmetric 2 x 2 (asym ignored)
  Eigen::MatrixXd m22(2, 2);
  m22 << 1, 2, 3, 4;
  stan::test::expect_ad(tols, f, m22);

  // challenging 2 x 2
  Eigen::MatrixXd a22(2, 2);
  a22 << 1, 2, 2, 1;
  tols.hessian_hessian_ = 5e-1;
  tols.hessian_fvar_hessian_ = 5e-1;
  tols.grad_hessian_grad_hessian_ = relative_tolerance(1e-1, 5e-1);
  stan::test::expect_ad(tols, f, a22);
}

TEST(MathMixMatFun, eigenvaluesSymLogDet) {
  // logdet(A) can be calculated using eigenvalues of matrix A
  // the derivative of logdet(A) should be inverse(A)
  // See stan-dev/math/issues/1803

  Eigen::MatrixXd a(4, 4);
  // Random symmetric matrix
  a << 1.8904,  0.7204, -0.1599,  1.2028,
       0.7204,  7.3394,  2.0895, -0.6151,
      -0.1599,  2.0895,  2.4601, -1.7219,
       1.2028, -0.6151, -1.7219,  4.5260;

  stan::math::matrix_d a_inv = stan::math::inverse(a);

  stan::math::matrix_v a_v(a);
  auto w = eigenvalues_sym(a_v);
  auto logdet = stan::math::sum(stan::math::log(w));

  stan::math::set_zero_all_adjoints();
  logdet.grad();

  ASSERT_TRUE(a_inv.val().isApprox(a_v.adj()));
}
