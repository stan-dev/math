#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, eltMultiply) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::elt_multiply(x, y);
  };

  Eigen::VectorXd a2(2);
  a2 << 2, 5;
  Eigen::VectorXd b2(2);
  b2 << 10, 100;
  stan::test::expect_ad(f, a2, b2);

  Eigen::RowVectorXd c2(2);
  c2 << 2, 5;
  Eigen::RowVectorXd d2(2);
  d2 << 10, 100;
  stan::test::expect_ad(f, c2, d2);

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  Eigen::MatrixXd e23(2, 3);
  e23 << 2, 5, 7, 13, 29, 112;
  Eigen::MatrixXd f23(2, 3);
  f23 << 1e1, 1e2, 1e3, 1e4, 1e5, 1e6;
  stan::test::expect_ad(tols, f, e23, f23);

  Eigen::VectorXd v0(0);
  Eigen::RowVectorXd rv0(0);
  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, v0, v0);
  stan::test::expect_ad(f, rv0, rv0);
  stan::test::expect_ad(f, m00, m00);
}
