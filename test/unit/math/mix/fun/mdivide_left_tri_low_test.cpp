#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(MathMixMatFun, mdivideLeftTriLow) {
  auto f = [](const auto& x, const auto& y) {
    // no need to symmetrize because only uses view
    return stan::math::mdivide_left_tri_low(x, y);
  };

  // size zero inputs
  Eigen::MatrixXd m00(0, 0);
  Eigen::MatrixXd m02(0, 2);
  Eigen::VectorXd v0(0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad(f, m00, m02);
  stan::test::expect_ad(f, m00, v0);

  // signature 1 of 2: matrix-matrix
  Eigen::MatrixXd aa(1, 1);
  aa << 1;
  Eigen::MatrixXd bb(1, 1);
  bb << 2;
  stan::test::expect_ad(f, aa, bb);

  // signature 2 of 2: matrix-vector
  Eigen::VectorXd cc(1);
  cc << 3;
  stan::test::expect_ad(f, aa, cc);

  Eigen::MatrixXd a(2, 2);
  a << 1, std::numeric_limits<double>::quiet_NaN(), -3, 5;
  Eigen::MatrixXd b(2, 2);
  b << 2, 5, 12, 109;
  Eigen::VectorXd c(2);
  c << 3, 5;
  stan::test::expect_ad(f, a, b);
  stan::test::expect_ad(f, a, c);

  Eigen::MatrixXd a2(2, 2);
  a2 << 1, 0, -3, 5;
  Eigen::MatrixXd b2(2, 2);
  b2 << 2, 5, 12, 109;
  stan::test::expect_ad(f, a2, b2);

  Eigen::MatrixXd y(3, 3);
  y << 1, 0, 0, 2, 3, 0, 4, 5, 6;
  Eigen::VectorXd z(3);
  z << 1, 2, 3;
  Eigen::MatrixXd u(3, 3);
  u << 1, 2, 3, 6, 5, 4, 7, 8, 9;
  stan::test::expect_ad(f, y, z);
  stan::test::expect_ad(f, u, y);

  // exception cases
  Eigen::MatrixXd d(3, 2);
  d << 1, 2, 5, -19, 73, 31;
  Eigen::VectorXd e(3);
  e << 1, 2, 9;
  stan::test::expect_ad(f, d, b);
  stan::test::expect_ad(f, d, c);
  stan::test::expect_ad(f, a, e);

  Eigen::MatrixXd m33 = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd m44 = Eigen::MatrixXd::Zero(4, 4);
  Eigen::VectorXd v4 = Eigen::VectorXd::Zero(4);
  Eigen::RowVectorXd rv3 = Eigen::RowVectorXd::Zero(3);

  // exceptions: wrong sizes
  stan::test::expect_ad(f, m33, m44);
  stan::test::expect_ad(f, m33, v4);

  // exceptions: wrong types
  stan::test::expect_ad(f, m33, rv3);
}
