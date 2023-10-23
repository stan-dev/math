#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, eigendecomposeSym) {
  auto g = [](const auto& y) {
    // eigenvectors
    // maintain symmetry for finite diffs; ignore if not square
    if (y.rows() != y.cols()) {
      return std::get<0>(stan::math::eigendecompose_sym(y));
    }
    auto a = ((y + y.transpose()) * 0.5).eval();
    return std::get<0>(stan::math::eigendecompose_sym(a));
  };

  auto f = [](const auto& y) {
    // eigenvalues
    if (y.rows() != y.cols()) {
      return std::get<1>(stan::math::eigendecompose_sym(y));
    }
    auto a = ((y + y.transpose()) * 0.5).eval();
    return std::get<1>(stan::math::eigendecompose_sym(a));
  };

  stan::test::ad_tolerances values_tols;
  values_tols.hessian_hessian_ = 5e-2;
  values_tols.hessian_fvar_hessian_ = 5e-2;

  // 1 x 1
  Eigen::MatrixXd a11(1, 1);
  a11 << 2.3;
  stan::test::expect_ad(values_tols, f, a11);
  stan::test::expect_ad(g, a11);

  Eigen::MatrixXd a23(2, 3);
  a23 << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(values_tols, f, a23);
  stan::test::expect_ad(g, a23);

  Eigen::MatrixXd m22(2, 2);
  m22 << 1, 2, 3, 4;
  stan::test::expect_ad(values_tols, f, m22);
  stan::test::expect_ad(g, m22);

  Eigen::MatrixXd a33(3, 3);
  a33 << 1, 2, 3, 2, 5, 7.9, 3, 7.9, 1.08;
  stan::test::expect_ad(values_tols, f, a33);
  stan::test::expect_ad(g, a33);
}

TEST(MathMixMatFun, eigendecomposeSym_varmat) {
  auto g = [](const auto& y) {
    // eigenvectors
    // maintain symmetry for finite diffs; ignore if not square
    if (y.rows() != y.cols()) {
      return std::get<0>(stan::math::eigendecompose_sym(y));
    }
    auto a = ((y + y.transpose()) * 0.5).eval();
    return std::get<0>(stan::math::eigendecompose_sym(a));
  };

  auto f = [](const auto& y) {
    // eigenvalues
    if (y.rows() != y.cols()) {
      return std::get<1>(stan::math::eigendecompose_sym(y));
    }
    auto a = ((y + y.transpose()) * 0.5).eval();
    return std::get<1>(stan::math::eigendecompose_sym(a));
  };

  stan::test::ad_tolerances values_tols;
  values_tols.hessian_hessian_ = 5e-2;
  values_tols.hessian_fvar_hessian_ = 5e-2;

  stan::test::ad_tolerances vec_tols;
  vec_tols.hessian_hessian_ = 5e-2;
  vec_tols.hessian_fvar_hessian_ = 5e-2;

  Eigen::MatrixXd a11(1, 1);
  a11 << 2.3;
  stan::test::expect_ad(values_tols, f, a11);
  stan::test::expect_ad_matvar(values_tols, f, a11);
  stan::test::expect_ad(g, a11);
  stan::test::expect_ad_matvar(vec_tols, g, a11);

  Eigen::MatrixXd m22(2, 2);
  m22 << 1, 2, 3, 4;
  stan::test::expect_ad(values_tols, f, m22);
  stan::test::expect_ad_matvar(values_tols, f, m22);
  stan::test::expect_ad(g, m22);
  stan::test::expect_ad_matvar(vec_tols, g, m22);

  Eigen::MatrixXd a23(2, 3);
  a23 << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(values_tols, f, a23);
  stan::test::expect_ad_matvar(values_tols, f, a23);
  stan::test::expect_ad(g, a23);
  stan::test::expect_ad_matvar(vec_tols, g, a23);

  Eigen::MatrixXd a33(3, 3);
  a33 << 1, 2, 3, 2, 5, 7.9, 3, 7.9, 1.08;
  stan::test::expect_ad(values_tols, f, a33);
  stan::test::expect_ad_matvar(values_tols, f, a33);
  stan::test::expect_ad(g, a33);
  vec_tols.hessian_fvar_hessian_ = 1e-2;
  vec_tols.hessian_hessian_ = 1e-2;
  stan::test::expect_ad_matvar(vec_tols, g, a33);
}
