#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(MathMixMatFun, add) {
  auto f = [](const auto& x, const auto& y) { return stan::math::add(x, y); };

  double d = 2;
  Eigen::MatrixXd m11(1, 1);
  m11 << 1;
  Eigen::MatrixXd m11b(1, 1);
  m11b << 10;
  Eigen::VectorXd v1(1);
  v1 << 2;
  Eigen::VectorXd v1b(1);
  v1b << -20;
  Eigen::RowVectorXd rv1(1);
  rv1 << 3;
  Eigen::RowVectorXd rv1b(1);
  rv1b << -3;
  stan::test::expect_ad(f, d, v1);
  stan::test::expect_ad(f, d, rv1);
  stan::test::expect_ad(f, d, m11);
  stan::test::expect_ad(f, v1, d);
  stan::test::expect_ad(f, rv1, d);
  stan::test::expect_ad(f, m11, d);
  stan::test::expect_ad(f, v1, v1b);
  stan::test::expect_ad(f, rv1, rv1b);
  stan::test::expect_ad(f, m11, m11b);

  Eigen::MatrixXd m22(2, 2);
  m22 << 1, 2, 3, 4;
  Eigen::MatrixXd m22b(2, 2);
  m22b << 10, 20, 40, 100;
  Eigen::VectorXd v2(2);
  v2 << 10, 12;
  Eigen::VectorXd v2b(2);
  v2b << -1, -5;
  Eigen::RowVectorXd rv2(2);
  rv2 << 11, 15;
  Eigen::RowVectorXd rv2b(2);
  rv2b << 100, -1001;
  stan::test::expect_ad(f, d, v2);
  stan::test::expect_ad(f, d, rv2);
  stan::test::expect_ad(f, d, m22);
  stan::test::expect_ad(f, v2, d);
  stan::test::expect_ad(f, rv2, d);
  stan::test::expect_ad(f, m22, d);
  stan::test::expect_ad(f, v2, v2b);
  stan::test::expect_ad(f, rv2, rv2b);
  stan::test::expect_ad(f, m22b, m22b);

  Eigen::VectorXd v5(5);
  v5 << 1, 2, 3, 4, 5;
  Eigen::VectorXd v5b(5);
  v5 << 2, 3, 4, 5, 6;
  Eigen::VectorXd rv5(5);
  v5 << 1, 2, 3, 4, 5;
  Eigen::VectorXd rv5b(5);
  v5 << 2, 3, 4, 5, 6;
  stan::test::expect_ad(f, v5, v5b);
  stan::test::expect_ad(f, rv5, rv5b);

  Eigen::MatrixXd m22c(2, 2);
  m22c << -10, 1, 10, 0;
  Eigen::MatrixXd m22d(2, 2);
  m22d << 10, -10, 1, 2;
  stan::test::expect_ad(f, m22c, m22d);

  // these will throw
  stan::test::expect_ad(f, v1, v2);
  stan::test::expect_ad(f, rv1, rv2);
  stan::test::expect_ad(f, m11, m22);
}
