#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, eigenvectorsSym) {
  auto f = [](const auto& y) {
    // maintain symmetry for finite diffs; ignore if not square
    if (y.rows() != y.cols()) {
      return stan::math::eigenvectors_sym(y);
    }
    auto a = ((y + y.transpose()) * 0.5).eval();
    return stan::math::eigenvectors_sym(a);
  };

  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, m00);

  Eigen::MatrixXd a11(1, 1);
  a11 << 2;
  stan::test::expect_ad(f, a11);

  Eigen::MatrixXd m22(2, 2);
  m22 << 1, 2, 3, 4;
  stan::test::expect_ad(f, m22);

  Eigen::MatrixXd a23(2, 3);
  a23 << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(f, a23);

  stan::test::ad_tolerances tols;
  tols.hessian_fvar_hessian_ = 1e-2;
  tols.hessian_hessian_ = 1e-2;

  Eigen::MatrixXd a33(3, 3);
  a33 << 1, 2, 3, 2, 5, 7.9, 3, 7.9, 1.08;
  stan::test::expect_ad(tols, f, a33);
}
