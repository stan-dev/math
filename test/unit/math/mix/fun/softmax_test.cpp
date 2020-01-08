#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, softmax) {
  auto f = [](const auto& x) { return stan::math::softmax(x); };

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  Eigen::VectorXd a(0);
  stan::test::expect_ad(tols, f, a);

  Eigen::VectorXd b(1);
  b << 0;
  stan::test::expect_ad(tols, f, b);

  Eigen::VectorXd c(2);
  c << -1, 1;
  stan::test::expect_ad(tols, f, c);

  Eigen::VectorXd d(3);
  d << -1, 1, 10;
  stan::test::expect_ad(tols, f, d);

  Eigen::VectorXd d2(3);
  d2 << 0.5, -1, 3;
  stan::test::expect_ad(tols, f, d2);

  Eigen::VectorXd d3(3);
  d3 << 4, 3, -2;
  stan::test::expect_ad(tols, f, d3);

  Eigen::VectorXd d4(3);
  d4 << 0, 3, -1;
  stan::test::expect_ad(tols, f, d4);
}
