#include <test/unit/math/test_ad.hpp>

void expect_upper_sym_tri(const Eigen::MatrixXd& m) {
  auto f = [](const auto& x) { 
		return stan::math::symmetrize_from_upper_tri(x); 
	};
  Eigen::VectorXd v = stan::test::to_vector(m);
  Eigen::RowVectorXd rv = stan::test::to_row_vector(m);
  stan::test::expect_ad(f, m);
  stan::test::expect_ad(f, v);
  stan::test::expect_ad(f, rv);
}

TEST(MathMixMatFun, symmetrize_from_upper_tri) {
  Eigen::MatrixXd a(0, 0);
  expect_upper_sym_tri(a);

  Eigen::MatrixXd b(1, 1);
  b << -1.2;
  expect_upper_sym_tri(b);

  Eigen::MatrixXd c(2, 2);
  c << -1, 2, 5, 10;
  expect_upper_sym_tri(c);

  Eigen::MatrixXd d(2, 3);
  d << -1, 2, -3, 5, 10, 100;
  expect_upper_sym_tri(d);

  Eigen::MatrixXd e(3, 1);
  e << 1, 2, 3;
  expect_upper_sym_tri(e);

  Eigen::MatrixXd g(3, 3);
  g << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  expect_upper_sym_tri(g);
}