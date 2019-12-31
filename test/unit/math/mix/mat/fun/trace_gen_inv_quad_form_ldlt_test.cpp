#include <test/unit/math/test_ad.hpp>

template <typename T>
stan::math::LDLT_factor<T, -1, -1> to_ldlt(const Eigen::Matrix<T, -1, -1>& a) {
  if (a.size() == 0) {
    return {};
  }

  stan::math::LDLT_factor<T, -1, -1> ldlt_a;
  ldlt_a.compute(a);
  return ldlt_a;
}

TEST(mathMixMatFun, traceGenInvQuadForm) {
  auto f = [](const auto& c, const auto& a, const auto& b) {
    auto ldlt_a = to_ldlt(a);
    return stan::math::trace_gen_inv_quad_form_ldlt(c, ldlt_a, b);
  };

  Eigen::MatrixXd a00(0, 0);
  Eigen::MatrixXd b00(0, 0);
  Eigen::MatrixXd c00(0, 0);
  stan::test::expect_ad(f, c00, a00, b00);

  Eigen::MatrixXd b02(0, 2);
  stan::test::expect_ad(f, c00, a00, b02);

  Eigen::MatrixXd a11(1, 1);
  a11 << 1;
  Eigen::MatrixXd b11(1, 1);
  b11 << 2;
  Eigen::MatrixXd c11(1, 1);
  c11 << -3;
  stan::test::expect_ad(f, c11, a11, b11);

  // tolerance very low for gradients with var; tolerance OK for fvar
  // which uses templated prim implementation using templated prim
  // impl for rev still has issues below gradient_grad_ = 1e0; hessian
  // tolerance also low
  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 1e3;
  tols.hessian_hessian_ = 1e0;
  tols.hessian_fvar_hessian_ = 1e0;

  Eigen::MatrixXd a(4, 4);
  a << 2, 3, 4, 5, 6, 10, 2, 2, 7, 2, 7, 1, 8, 2, 1, 112;
  Eigen::MatrixXd b(4, 2);
  b << 100, 10, 0, 1, -3, -3, 5, 2;
  Eigen::MatrixXd c(2, 2);
  c.setIdentity(2, 2);
  stan::test::expect_ad(tols, f, c, a, b);

  Eigen::MatrixXd d(2, 2);
  d << 1, 2, 3, 4;
  stan::test::expect_ad(tols, f, d, a, b);

  Eigen::MatrixXd A(2, 2);
  A << 2, 3, 3, 7;
  Eigen::MatrixXd B(2, 2);
  B << 5, 6, 7, 8;
  Eigen::MatrixXd D(2, 2);
  D << 9, 10, 11, 12;
  stan::test::expect_ad(tols, f, D, A, B);

  // exception tests
  // non-square second arg
  Eigen::MatrixXd a34(3, 4);
  stan::test::expect_ad(f, c, a34, b);

  // non-square first arg
  Eigen::MatrixXd c23(2, 3);
  stan::test::expect_ad(f, c23, a, b);

  // a, b not multiplicable
  Eigen::MatrixXd b32(3, 2);
  stan::test::expect_ad(f, c, a, b32);

  // b, c not multiplicable
  Eigen::MatrixXd b43(4, 3);
  stan::test::expect_ad(f, c, a, b43);
}
