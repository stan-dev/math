#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(MathMixMatFun, mdivideRightTriLow) {
  auto f = [](const auto& x, const auto& y) {
    // no need to symmetrize because only uses view
    return stan::math::mdivide_right_tri_low(x, y);
  };

  // size zero inputs
  Eigen::MatrixXd m00(0, 0);
  Eigen::MatrixXd m20(2, 0);
  Eigen::RowVectorXd rv0(0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad(f, m20, m00);
  stan::test::expect_ad(f, rv0, m00);

  // signature 1 of 2: matrix / matrix
  Eigen::MatrixXd aa(1, 1);
  aa << 1;
  Eigen::MatrixXd bb(1, 1);
  bb << 2;
  stan::test::expect_ad(f, aa, bb);

  // signature 2 of 2: row-vector / matrix
  Eigen::RowVectorXd cc(1);
  cc << 3;
  stan::test::expect_ad(f, cc, aa);

  Eigen::MatrixXd a(2, 2);
  a << 1, std::numeric_limits<double>::quiet_NaN(), -3, 5;
  Eigen::MatrixXd b(2, 2);
  b << 2, 5, 12, 109;
  Eigen::RowVectorXd c(2);
  c << 3, 5;
  stan::test::expect_ad(f, b, a);
  stan::test::expect_ad(f, c, a);

  Eigen::MatrixXd a2(2, 2);
  a2 << 1, 0, -3, 5;
  Eigen::MatrixXd b2(2, 2);
  b2 << 2, 5, 12, 109;
  stan::test::expect_ad(f, b2, a2);

  Eigen::MatrixXd y(3, 3);
  y << 1, 0, 0, 2, 3, 0, 4, 5, 6;
  Eigen::RowVectorXd z(3);
  z << 1, 2, 3;
  Eigen::MatrixXd u(3, 3);
  u << 1, 2, 3, 6, 5, 4, 7, 8, 9;
  stan::test::expect_ad(f, z, y);
  stan::test::expect_ad(f, u, y);

  // exception cases
  Eigen::MatrixXd d(3, 2);
  d << 1, 2, 5, -19, 73, 31;
  Eigen::RowVectorXd e(3);
  e << 1, 2, 9;
  stan::test::expect_ad(f, b, d);
  stan::test::expect_ad(f, c, d);
  stan::test::expect_ad(f, e, a);

  Eigen::MatrixXd m33 = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd m44 = Eigen::MatrixXd::Zero(4, 4);
  Eigen::RowVectorXd rv4 = Eigen::VectorXd::Zero(4);
  Eigen::VectorXd v3 = Eigen::RowVectorXd::Zero(3);

  // exceptions: wrong sizes
  stan::test::expect_ad(f, m33, m44);
  stan::test::expect_ad(f, rv4, m33);

  // exceptions: wrong types
  stan::test::expect_ad(f, v3, m33);
}
