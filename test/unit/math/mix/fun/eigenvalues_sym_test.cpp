#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, eigenvaluesSym) {
  auto f = [](const auto& y) {
    // need to maintain symmetry for finite diffs
    auto a = ((y + y.transpose()) * 0.5).eval();
    return stan::math::eigenvalues_sym(a);
  };

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 5e-2;
  tols.hessian_fvar_hessian_ = 5e-2;

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
  tols.grad_hessian_grad_hessian_ = stan::test::relative_tolerance(1e-1, 5e-1);
  stan::test::expect_ad(tols, f, a22);
}

TEST(MathMixMatFun, eigenvaluesSym_varmat) {
  auto f = [](const auto& y) {
    // need to maintain symmetry for finite diffs
    auto a = stan::math::multiply((y + y.transpose()), 0.5).eval();
    return stan::math::eigenvalues_sym(a);
  };

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 5e-2;
  tols.hessian_fvar_hessian_ = 5e-2;

  Eigen::MatrixXd a11(1, 1);
  a11 << 2.3;
  stan::test::expect_ad(tols, f, a11);
  stan::test::expect_ad_matvar(tols, f, a11);

  Eigen::MatrixXd a33(3, 3);
  a33 << 1, 2, 3, 2, 5, 7.9, 3, 7.9, 1.08;
  stan::test::expect_ad(tols, f, a33);
  stan::test::expect_ad_matvar(tols, f, a33);

  Eigen::MatrixXd m22(2, 2);
  m22 << 1, 2, 3, 4;
  stan::test::expect_ad(tols, f, m22);
  stan::test::expect_ad_matvar(tols, f, m22);

  Eigen::MatrixXd a22(2, 2);
  a22 << 1, 2, 2, 1;
  tols.hessian_hessian_ = 5e-1;
  tols.hessian_fvar_hessian_ = 5e-1;
  tols.grad_hessian_grad_hessian_ = stan::test::relative_tolerance(1e-1, 5e-1);
  stan::test::expect_ad(tols, f, a22);
  stan::test::expect_ad_matvar(tols, f, a22);
}
