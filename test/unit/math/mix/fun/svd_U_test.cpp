#include <test/unit/math/test_ad.hpp>
//#include <stan/math/rev/core.hpp>
TEST(MathMixMatFun, svd_U) {
  auto f = [](const auto& x) { return stan::math::svd_U(x); };

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, m00);

  Eigen::MatrixXd m11(1, 1);
  m11 << 1.1;
  stan::test::expect_ad(tols, f, m11);

  Eigen::MatrixXd m22(2, 2);
  m22 << 3, -5, 7, 11;
  stan::test::expect_ad(tols, f, m22);

  Eigen::MatrixXd m23(2, 3);
  m23 << 3, 5, -7, -11, 13, -17;
  stan::test::expect_ad(tols, f, m23);

  Eigen::MatrixXd m32(3, 2);
  m32 << 1, 3, -5, 7, 9, -11;
  stan::test::expect_ad(tols, f, m32);

  Eigen::MatrixXd a22(2, 2);
  a22 << 1, 2, 3, 4;
  stan::test::expect_ad(tols, f, a22);
}
