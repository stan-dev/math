#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, appendRow) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::append_row(x, y);
  };

  Eigen::MatrixXd m22(2, 2);
  m22 << 2, 3, 9, -1;
  Eigen::MatrixXd m22b(2, 2);
  m22b << 4, 3, 0, 1;
  Eigen::MatrixXd m32(3, 2);
  m32 << 4, 3, 0, 1, 11, 12;
  Eigen::VectorXd v2(2);
  v2 << 17, -32;
  Eigen::VectorXd v2b(2);
  v2b << 11, 9;
  Eigen::VectorXd rv2(2);
  rv2 << 17, -32;
  Eigen::VectorXd rv2b(2);
  rv2b << 11, 9;
  Eigen::RowVectorXd v4(4);
  v4 << 2, 3, 9, -1;
  Eigen::RowVectorXd v3(3);
  v3 << 4, 3, 0.4;

  // (matrix, matrix)
  stan::test::expect_ad(f, m22, m22b);
  stan::test::expect_ad(f, m32, m22);
  stan::test::expect_ad_matvar(f, m22, m22b);
  stan::test::expect_ad_matvar(f, m32, m22);

  // (row_vector, matrix)
  stan::test::expect_ad(f, rv2, m22);
  stan::test::expect_ad_matvar(f, rv2, m22);

  // (matrix, row_vector)
  stan::test::expect_ad(f, m22, rv2);
  stan::test::expect_ad_matvar(f, m22, rv2);

  // (vector, vector)
  stan::test::expect_ad(f, v4, v3);
  stan::test::expect_ad(f, v3, v4);
  stan::test::expect_ad_matvar(f, v4, v3);
  stan::test::expect_ad_matvar(f, v3, v4);

  // (row vector, row vector)
  stan::test::expect_ad(f, rv2, rv2b);
  stan::test::expect_ad_matvar(f, rv2, rv2b);

  Eigen::MatrixXd m00(0, 0);
  Eigen::MatrixXd m30(3, 0);
  Eigen::VectorXd v0(0);
  Eigen::RowVectorXd rv0(0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad(f, m30, m00);
  stan::test::expect_ad(f, m00, rv0);
  stan::test::expect_ad(f, v0, v0);
  stan::test::expect_ad(f, rv0, rv0);

  stan::test::expect_ad_matvar(f, m00, m00);
  stan::test::expect_ad_matvar(f, m30, m00);
  stan::test::expect_ad_matvar(f, m00, rv0);
  stan::test::expect_ad_matvar(f, v0, v0);
  stan::test::expect_ad_matvar(f, rv0, rv0);
}
