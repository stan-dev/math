#include <test/unit/math/test_ad.hpp>
#include <test/unit/util.hpp>

TEST(MathMixMatFun, softmax) {
  using stan::test::expect_ad_matvar;
  auto f = [](const auto& x) { return stan::math::softmax(x); };

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-2;
  tols.hessian_fvar_hessian_ = 1e-2;

  // Column vectors
  Eigen::VectorXd a(0);
  stan::test::expect_ad(tols, f, a);
  expect_ad_matvar(f, a);
  Eigen::VectorXd b(1);
  b << 0;
  stan::test::expect_ad(tols, f, b);
  expect_ad_matvar(f, b);

  Eigen::VectorXd c(2);
  c << -1, 1;
  stan::test::expect_ad(tols, f, c);
  expect_ad_matvar(f, c);

  Eigen::VectorXd d(3);
  d << -1, 1, 10;
  stan::test::expect_ad(tols, f, d);
  expect_ad_matvar(f, d);

  Eigen::VectorXd d2(3);
  d2 << 0.5, -1, 3;
  stan::test::expect_ad(tols, f, d2);
  expect_ad_matvar(f, d2);

  Eigen::VectorXd d3(3);
  d3 << 4, 3, -2;
  stan::test::expect_ad(tols, f, d3);
  expect_ad_matvar(f, d3);

  Eigen::VectorXd d4(3);
  d4 << 0, 3, -1;
  stan::test::expect_ad(tols, f, d4);
  expect_ad_matvar(f, d4);

  // Row vectors
  Eigen::RowVectorXd ra(0);
  stan::test::expect_ad(tols, f, ra);
  expect_ad_matvar(f, ra);

  Eigen::RowVectorXd rb(1);
  rb << 0;
  stan::test::expect_ad(tols, f, rb);
  expect_ad_matvar(f, rb);

  Eigen::RowVectorXd rc(2);
  rc << -1, 1;
  stan::test::expect_ad(tols, f, rc);
  expect_ad_matvar(f, rc);

  Eigen::RowVectorXd rd(3);
  rd << -1, 1, 10;
  stan::test::expect_ad(tols, f, rd);
  expect_ad_matvar(f, rd);

  Eigen::RowVectorXd rd2(3);
  rd2 << 0.5, -1, 3;
  stan::test::expect_ad(tols, f, rd2);
  expect_ad_matvar(f, rd2);
}
