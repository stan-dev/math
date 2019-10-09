#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixMatFun, inverse) {
  auto f = [](const auto& x) { return stan::math::inverse(x); };

  Eigen::MatrixXd u(1, 1);
  u << 2;
  stan::test::expect_ad(f, u);

  Eigen::MatrixXd v(2, 2);
  v << 2, 3, 5, 7;
  stan::test::expect_ad(f, v);

  // issues around zero require looser tolerances for hessians
  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 0.01;
  tols.hessian_fvar_hessian_ = 0.01;

  Eigen::MatrixXd w(3, 3);
  w << 2, 3, 5, 7, 11, 13, 17, 19, 23;
  stan::test::expect_ad(tols, f, w);

  // even worse here in some cases
  stan::test::ad_tolerances tols2;
  tols2.hessian_hessian_ = 3.0;
  tols2.hessian_fvar_hessian_ = 3.0;
  Eigen::MatrixXd x(4, 4);
  x << 2, 3, 4, 5, 9, -1, 2, 2, 4, 3, 7, -1, 0, 1, 19, 112;
  stan::test::expect_ad(tols2, f, x);
}
