#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, subtract_transpose_test) {
  auto f = [](const auto& x) { return x - x.transpose(); };

  Eigen::MatrixXd x(2, 2);

  stan::test::expect_ad_matvar(f, x);
}

TEST(MathMixMatFun, subtract_block_test) {
  auto f = [](const auto& y, const auto& x) {
    return stan::math::subtract(y, x.block(0, 0, 2, 2));
  };

  Eigen::MatrixXd x(3, 3);

  stan::test::expect_ad_matvar(f, 1.0, x);
}

TEST(MathMixMatFun, subtract_1) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::subtract(x, y); };

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
}

TEST(MathMixMatFun, subtract_scal) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::subtract(x, y); };

  stan::test::expect_ad(f, 2, 5);
  stan::test::expect_ad(f, 2.0, 5);
  stan::test::expect_ad(f, 2, 5.0);
  stan::test::expect_ad(f, 2.0, 5.0);
}

TEST(MathMixMatFun, subtract_empty) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::subtract(x, y); };

  double d = 2;
  Eigen::MatrixXd m00(0, 0);
  Eigen::MatrixXd m00b(0, 0);
  Eigen::VectorXd v0(0);
  Eigen::VectorXd v0b(0);
  Eigen::RowVectorXd rv0(0);
  Eigen::RowVectorXd rv0b(0);
  stan::test::expect_ad(f, d, v0);
  stan::test::expect_ad(f, d, rv0);
  stan::test::expect_ad(f, d, m00);
  stan::test::expect_ad(f, v0, d);
  stan::test::expect_ad(f, rv0, d);
  stan::test::expect_ad(f, m00, d);
  stan::test::expect_ad(f, v0, v0b);
  stan::test::expect_ad(f, rv0, rv0b);
  stan::test::expect_ad(f, m00, m00b);
}

TEST(MathMixMatFun, subtract_scalar_mat) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::subtract(x, y); };

  double d = 2;
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
}

TEST(MathMixMatFun, subtract_vec_mat) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::subtract(x, y); };

  Eigen::VectorXd v5(5);
  v5 << 1, 2, 3, 4, 5;
  Eigen::VectorXd v5b(5);
  v5b << 2, 3, 4, 5, 6;
  Eigen::VectorXd rv5(5);
  rv5 << 1, 2, 3, 4, 5;
  Eigen::VectorXd rv5b(5);
  rv5b << 2, 3, 4, 5, 6;
  stan::test::expect_ad(f, v5, v5b);
  stan::test::expect_ad(f, rv5, rv5b);
}

TEST(MathMixMatFun, subtract_mat_mat) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::subtract(x, y); };

  Eigen::MatrixXd m22c(2, 2);
  m22c << -10, 1, 10, 0;
  Eigen::MatrixXd m22d(2, 2);
  m22d << 10, -10, 1, 2;
  stan::test::expect_ad(f, m22c, m22d);

  Eigen::MatrixXd m23(2, 3);
  m23 << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd m23b(2, 3);
  m23b << 10, 12, 14, -1, -3, -10;
  stan::test::expect_ad(f, m23, m23b);
}

// these will throw
TEST(MathMixMatFun, subtract_throw) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::subtract(x, y); };

  Eigen::VectorXd v1(1);
  v1 << 2;
  Eigen::VectorXd v2(2);
  v2 << 10, 12;
  stan::test::expect_ad(f, v1, v2);

  Eigen::RowVectorXd rv1(1);
  rv1 << 3;
  Eigen::RowVectorXd rv2(2);
  rv2 << 11, 15;
  stan::test::expect_ad(f, rv1, rv2);

  Eigen::MatrixXd m22(2, 2);
  m22 << 1, 2, 3, 4;
  Eigen::MatrixXd m11(1, 1);
  m11 << 1;
  stan::test::expect_ad(f, m11, m22);
}
