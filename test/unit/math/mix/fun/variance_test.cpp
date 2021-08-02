#include <test/unit/math/test_ad.hpp>
#include <vector>

void expect_variance(const Eigen::MatrixXd& m) {
  auto f = [](const auto& x) { return stan::math::variance(x); };
  Eigen::VectorXd v = stan::test::to_vector(m);
  Eigen::RowVectorXd rv = stan::math::to_row_vector(m);
  std::vector<double> sv = stan::math::to_array_1d(m);

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = 1e-1;
  tols.hessian_fvar_hessian_ = 1e-1;
  stan::test::expect_ad(tols, f, m);
  stan::test::expect_ad(tols, f, v);
  stan::test::expect_ad(tols, f, rv);
  stan::test::expect_ad(tols, f, sv);

  stan::test::expect_ad_matvar(tols, f, m);
  stan::test::expect_ad_matvar(tols, f, v);
  stan::test::expect_ad_matvar(tols, f, rv);
}

TEST(MathMixMatFun, variance) {
  Eigen::MatrixXd a(0, 0);
  expect_variance(a);

  Eigen::MatrixXd b(1, 1);
  b << -1.2;
  expect_variance(b);

  Eigen::MatrixXd b11(1, 1);
  b11 << 12.9;
  expect_variance(b11);

  Eigen::MatrixXd c(2, 2);
  c << -1, 2, 5, 10;
  expect_variance(c);

  Eigen::MatrixXd c21(2, 1);
  c21 << 3, 1.7;
  expect_variance(c21);

  Eigen::MatrixXd d(2, 3);
  d << -1, 2, -3, 5, 10, 100;
  expect_variance(d);

  Eigen::MatrixXd d23(2, 3);
  d << 1, 2, 3, 4, 5, 6;
  expect_variance(d);

  Eigen::MatrixXd e(3, 1);
  e << 1, 2, 3;
  expect_variance(e);

  Eigen::MatrixXd e31(3, 1);
  e31 << 0.5, 2, 3.5;
  expect_variance(e31);

  Eigen::MatrixXd g(3, 3);
  g << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  expect_variance(g);

  Eigen::MatrixXd h(4, 1);
  h << 1, 2, 5, 5;
  expect_variance(h);
}
