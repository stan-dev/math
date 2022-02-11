#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, traceGenQuadForm) {
  using stan::test::relative_tolerance;
  auto f = [](const auto& x1, const auto& x2, const auto& x3) {
    return stan::math::trace_gen_quad_form(x1, x2, x3);
  };

  Eigen::MatrixXd a00(0, 0);
  Eigen::MatrixXd b00(0, 0);
  Eigen::MatrixXd c00(0, 0);
  stan::test::expect_ad(f, c00, a00, b00);
  stan::test::expect_ad_matvar(f, c00, a00, b00);

  Eigen::MatrixXd a11(1, 1);
  a11 << 1;
  Eigen::MatrixXd b11(1, 1);
  b11 << 2;
  Eigen::MatrixXd c11(1, 1);
  c11 << -3;
  stan::test::expect_ad(f, c11, a11, b11);
  stan::test::expect_ad_matvar(f, c11, a11, b11);

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = relative_tolerance(0.1, 0.15);
  tols.hessian_fvar_hessian_ = relative_tolerance(0.1, 0.15);

  Eigen::MatrixXd a(4, 4);
  a << 2, 3, 4, 5, 6, 10, 2, 2, 7, 2, 7, 1, 8, 2, 1, 112;
  Eigen::MatrixXd b(4, 2);
  b << 100, 10, 0, 1, -3, -3, 5, 2;
  Eigen::MatrixXd c(2, 2);
  c.setIdentity(2, 2);
  stan::test::expect_ad(tols, f, c, a, b);
  stan::test::expect_ad_matvar(f, c, a, b);

  // exception tests
  // non-square second arg
  Eigen::MatrixXd a34(3, 4);
  stan::test::expect_ad(tols, f, c, a34, b);
  stan::test::expect_ad_matvar(tols, f, c, a34, b);

  // non-square first arg
  Eigen::MatrixXd c23(2, 3);
  stan::test::expect_ad(tols, f, c23, a, b);
  stan::test::expect_ad_matvar(tols, f, c23, a, b);

  // a, b not multiplicable
  Eigen::MatrixXd b32(3, 2);
  stan::test::expect_ad(tols, f, c, a, b32);
  stan::test::expect_ad_matvar(tols, f, c, a, b32);

  // b, c not multiplicable
  Eigen::MatrixXd b43(4, 3);
  stan::test::expect_ad(tols, f, c, a, b43);
  stan::test::expect_ad_matvar(tols, f, c, a, b43);
}
