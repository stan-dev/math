#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, traceGenInvQuadForm) {
  auto f = [](const auto& c, const auto& a, const auto& b) {
    auto x_sym = stan::math::multiply(0.5, a + a.transpose());
    auto ldlt_a = stan::math::make_ldlt_factor(x_sym);
    return stan::math::trace_gen_inv_quad_form_ldlt(c, ldlt_a, b);
  };

  Eigen::MatrixXd a00(0, 0);
  Eigen::MatrixXd b00(0, 0);
  Eigen::MatrixXd c00(0, 0);
  stan::test::expect_ad(f, c00, a00, b00);
  stan::test::expect_ad_matvar(f, c00, a00, b00);

  Eigen::MatrixXd b02(0, 2);
  stan::test::expect_ad(f, c00, a00, b02);
  stan::test::expect_ad_matvar(f, c00, a00, b02);

  Eigen::MatrixXd a11(1, 1);
  a11 << 1;
  Eigen::MatrixXd b11(1, 1);
  b11 << 2;
  Eigen::MatrixXd c11(1, 1);
  c11 << -3;
  stan::test::expect_ad(f, c11, a11, b11);
  stan::test::expect_ad_matvar(f, c11, a11, b11);

  // tolerance very low for gradients
  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = stan::test::relative_tolerance(0.1, 0.15);
  tols.hessian_hessian_ = stan::test::relative_tolerance(0.1, 0.15);
  tols.hessian_fvar_hessian_ = stan::test::relative_tolerance(0.1, 0.15);

  Eigen::MatrixXd a(4, 4);
  a << 6, 3, 4, 1, 3, 7, 2, 2, 4, 2, 7, 1, 1, 2, 1, 5;
  Eigen::MatrixXd b(4, 2);
  b << 100, 10, 0, 1, -3, -3, 5, 2;
  Eigen::MatrixXd c(2, 2);
  c.setIdentity(2, 2);
  stan::test::expect_ad(tols, f, c, a, b);
  stan::test::expect_ad_matvar(tols, f, c, a, b);

  Eigen::MatrixXd d(2, 2);
  d << 1, 2, 3, 4;
  stan::test::expect_ad(tols, f, d, a, b);
  stan::test::expect_ad_matvar(tols, f, d, a, b);

  Eigen::MatrixXd A(2, 2);
  A << 3, 1, 1, 4;
  Eigen::MatrixXd B(2, 2);
  B << 5, 6, 7, 8;
  Eigen::MatrixXd D(2, 2);
  D << 9, 10, 11, 12;
  stan::test::expect_ad(tols, f, D, A, B);
  stan::test::expect_ad_matvar(tols, f, D, A, B);

  // exception tests
  // non-square first arg
  Eigen::MatrixXd c23(2, 3);
  stan::test::expect_ad(f, c23, a, b);
  stan::test::expect_ad_matvar(f, c23, a, b);

  // a, b not multiplicable
  Eigen::MatrixXd b32(3, 2);
  stan::test::expect_ad(f, c, a, b32);
  stan::test::expect_ad_matvar(f, c, a, b32);

  // b, c not multiplicable
  Eigen::MatrixXd b43(4, 3);
  stan::test::expect_ad(f, c, a, b43);
  stan::test::expect_ad_matvar(f, c, a, b43);

  stan::math::recover_memory();
}

TEST(mathMixMatFun, traceGenInvQuadForm_vec) {
  auto f = [](const auto& c, const auto& a, const auto& b) {
    auto x_sym = stan::math::multiply(0.5, a + a.transpose());
    auto ldlt_a = stan::math::make_ldlt_factor(x_sym);
    return stan::math::trace_gen_inv_quad_form_ldlt(c, ldlt_a, b);
  };

  Eigen::MatrixXd a00(0, 0);
  Eigen::MatrixXd b00(0, 0);
  Eigen::VectorXd c0(0);
  stan::test::expect_ad(f, c0, a00, b00);
  stan::test::expect_ad_matvar(f, c0, a00, b00);

  Eigen::MatrixXd b02(0, 2);
  stan::test::expect_ad(f, c0, a00, b02);
  stan::test::expect_ad_matvar(f, c0, a00, b02);

  Eigen::MatrixXd a11(1, 1);
  a11 << 1;
  Eigen::MatrixXd b11(1, 1);
  b11 << 2;
  Eigen::VectorXd c1(1);
  c1 << -3;
  stan::test::expect_ad(f, c1, a11, b11);
  stan::test::expect_ad_matvar(f, c1, a11, b11);

  // tolerance very low for gradients
  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = stan::test::relative_tolerance(0.1, 0.15);
  tols.hessian_hessian_ = stan::test::relative_tolerance(0.1, 0.15);
  tols.hessian_fvar_hessian_ = stan::test::relative_tolerance(0.1, 0.15);

  Eigen::MatrixXd a(4, 4);
  a << 6, 3, 4, 1, 3, 7, 2, 2, 4, 2, 7, 1, 1, 2, 1, 5;
  Eigen::MatrixXd b(4, 2);
  b << 100, 10, 0, 1, -3, -3, 5, 2;
  Eigen::VectorXd c(2);
  c << 2.0, 0.5;
  stan::test::expect_ad(tols, f, c, a, b);
  stan::test::expect_ad_matvar(tols, f, c, a, b);

  Eigen::VectorXd d(2);
  d << 1, 3;
  stan::test::expect_ad(tols, f, d, a, b);
  stan::test::expect_ad_matvar(tols, f, d, a, b);

  Eigen::MatrixXd A(2, 2);
  A << 3, 1, 1, 4;
  Eigen::MatrixXd B(2, 2);
  B << 5, 6, 7, 8;
  Eigen::VectorXd D(2);
  D << 9, 10;
  stan::test::expect_ad(tols, f, D, A, B);
  stan::test::expect_ad_matvar(tols, f, D, A, B);

  // exception tests
  // a, b not multiplicable
  Eigen::MatrixXd b32(3, 2);
  stan::test::expect_ad(f, c, a, b32);
  stan::test::expect_ad_matvar(f, c, a, b32);

  // b, c not multiplicable
  Eigen::MatrixXd b43(4, 3);
  stan::test::expect_ad(f, c, a, b43);
  stan::test::expect_ad_matvar(f, c, a, b43);

  stan::math::recover_memory();
}
