#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, mdivideRightTri) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::mdivide_right_tri<Eigen::Lower>(x, y);
  };
  auto f_up = [](const auto& x, const auto& y) {
    return stan::math::mdivide_right_tri<Eigen::Upper>(x, y);
  };

  // size zero inputs
  Eigen::MatrixXd m00(0, 0);
  Eigen::RowVectorXd rv0(0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad(f, rv0, m00);

  // signature 1 of 2: matrix / matrix
  Eigen::MatrixXd aa(1, 1);
  aa << 1;
  Eigen::MatrixXd bb(1, 1);
  bb << 2;
  stan::test::expect_ad(f, aa, bb);
  stan::test::expect_ad(f_up, aa, bb);

  // signature 2 of 2: row-vect / matrix
  Eigen::RowVectorXd cc(1);
  cc << 3;
  stan::test::expect_ad(f, cc, aa);
  stan::test::expect_ad(f_up, cc, aa);

  Eigen::MatrixXd a(2, 2);
  a << 2, 0, 5, 7;
  Eigen::MatrixXd a_tr = a.transpose();
  Eigen::MatrixXd b(2, 2);
  b << 12, 13, 15, 17;
  Eigen::RowVectorXd c(2);
  c << 3, 5;
  stan::test::expect_ad(f, a, a);
  stan::test::expect_ad(f, b, a);
  stan::test::expect_ad(f, c, a);
  stan::test::expect_ad(f_up, a, a_tr);
  stan::test::expect_ad(f_up, b, a_tr);
  stan::test::expect_ad(f_up, c, a_tr);

  Eigen::MatrixXd y(3, 3);
  y << 1, 0, 0, 2, 3, 0, 4, 5, 6;
  Eigen::MatrixXd y_tr = y.transpose();
  Eigen::RowVectorXd z(3);
  z << 1, 2, 3;
  Eigen::MatrixXd u(3, 3);
  u << 1, 2, 3, 6, 5, 4, 7, 8, 9;
  stan::test::expect_ad(f, z, y);
  stan::test::expect_ad(f, y, u);
  stan::test::expect_ad(f_up, z, y_tr);
  stan::test::expect_ad(f_up, y, y_tr);

  Eigen::MatrixXd uu(2, 2);
  uu << 3, 0, 1, 4;
  Eigen::MatrixXd vv(2, 2);
  vv << 2, 3, 5, 7;
  stan::test::expect_ad(f, vv, uu);

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
  Eigen::RowVectorXd v3 = Eigen::VectorXd::Zero(3);
  Eigen::VectorXd rv4 = Eigen::RowVectorXd::Zero(4);

  // exceptions: wrong sizes
  stan::test::expect_ad(f, m33, m44);
  stan::test::expect_ad(f, rv4, m33);

  // exceptions: wrong types
  stan::test::expect_ad(f, v3, m33);
}
