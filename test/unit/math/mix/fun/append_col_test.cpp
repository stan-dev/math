#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, appendCol) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::append_col(x, y);
  };

  Eigen::MatrixXd m22(2, 2);
  m22 << 2, 3, 9, -1;
  Eigen::MatrixXd m22b(2, 2);
  m22b << 4, 3, 0, 1;
  Eigen::MatrixXd m23b(2, 3);
  m23b << 4, 3, 0, 1, 11, 12;
  Eigen::VectorXd v2(2);
  v2 << 17, -32;
  Eigen::VectorXd v2b(2);
  v2b << 11, 9;
  Eigen::RowVectorXd rv4(4);
  rv4 << 2, 3, 9, -1;
  Eigen::RowVectorXd rv3(3);
  rv3 << 4, 3, 0.4;

  // (matrix, matrix)
  stan::test::expect_ad(f, m22, m22b);
  stan::test::expect_ad(f, m22, m23b);
  stan::test::expect_ad_matvar(f, m22, m22b);
  stan::test::expect_ad_matvar(f, m22, m23b);

  // (vector, matrix)
  stan::test::expect_ad(f, v2, m22);
  stan::test::expect_ad_matvar(f, v2, m22);

  // (matrix, vector)
  stan::test::expect_ad(f, m22, v2);
  stan::test::expect_ad_matvar(f, m22, v2);

  // (vector, vector)
  stan::test::expect_ad(f, v2, v2b);
  stan::test::expect_ad_matvar(f, v2, v2b);

  // (row vector, row vector)
  stan::test::expect_ad(f, rv3, rv4);
  stan::test::expect_ad(f, rv4, rv3);
  stan::test::expect_ad_matvar(f, rv3, rv4);
  stan::test::expect_ad_matvar(f, rv4, rv3);

  Eigen::MatrixXd m00(0, 0);
  Eigen::MatrixXd m03(0, 3);
  Eigen::VectorXd v0(0);
  Eigen::RowVectorXd rv0(0);

  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad(f, m00, m03);
  stan::test::expect_ad(f, v0, m00);
  stan::test::expect_ad(f, v0, v0);
  stan::test::expect_ad(f, rv0, rv0);

  stan::test::expect_ad_matvar(f, m00, m00);
  stan::test::expect_ad_matvar(f, m00, m03);
  stan::test::expect_ad_matvar(f, v0, m00);
  stan::test::expect_ad_matvar(f, v0, v0);
  stan::test::expect_ad_matvar(f, rv0, rv0);
}
