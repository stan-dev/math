#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, mdivideLeftTri) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::mdivide_left_tri<Eigen::Lower>(x, y);
  };
  auto f_up = [](const auto& x, const auto& y) {
    return stan::math::mdivide_left_tri<Eigen::Upper>(x, y);
  };

  // size zero inputs
  Eigen::MatrixXd m00(0, 0);
  Eigen::VectorXd v0(0);
  stan::test::expect_ad(f, m00, v0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad_matvar(f, m00, v0);
  stan::test::expect_ad_matvar(f, m00, m00);

  // signature 1 of 2: matrix-matrix
  Eigen::MatrixXd aa(1, 1);
  aa << 1;
  Eigen::MatrixXd bb(1, 1);
  bb << 2;
  stan::test::expect_ad(f, aa, bb);
  stan::test::expect_ad(f_up, aa, bb);
  stan::test::expect_ad_matvar(f, aa, bb);
  stan::test::expect_ad_matvar(f_up, aa, bb);

  // signature 2 of 2: matrix-vector
  Eigen::VectorXd cc(1);
  cc << 3;
  stan::test::expect_ad(f, aa, cc);
  stan::test::expect_ad(f_up, aa, cc);
  stan::test::expect_ad_matvar(f, aa, cc);
  stan::test::expect_ad_matvar(f_up, aa, cc);

  Eigen::MatrixXd a(2, 2);
  a << 2, 0, 5, 7;
  Eigen::MatrixXd a_tr = a.transpose();
  Eigen::MatrixXd b(2, 2);
  b << 12, 13, 15, 17;
  Eigen::VectorXd c(2);
  c << 3, 5;
  stan::test::expect_ad(f, a, a);
  stan::test::expect_ad(f, a, b);
  stan::test::expect_ad(f, a, c);
  stan::test::expect_ad(f_up, a_tr, a);
  stan::test::expect_ad(f_up, a_tr, b);
  stan::test::expect_ad(f_up, a_tr, c);
  stan::test::expect_ad_matvar(f, a, a);
  stan::test::expect_ad_matvar(f, a, b);
  stan::test::expect_ad_matvar(f, a, c);
  stan::test::expect_ad_matvar(f_up, a_tr, a);
  stan::test::expect_ad_matvar(f_up, a_tr, b);
  stan::test::expect_ad_matvar(f_up, a_tr, c);

  Eigen::MatrixXd y(3, 3);
  y << 1, 0, 0, 2, 3, 0, 4, 5, 6;
  Eigen::MatrixXd y_tr = y.transpose();
  Eigen::VectorXd z(3);
  z << 1, 2, 3;
  Eigen::MatrixXd u(3, 3);
  u << 1, 2, 3, 6, 5, 4, 7, 8, 9;
  stan::test::expect_ad(f, y, z);
  stan::test::expect_ad(f, u, y);
  stan::test::expect_ad(f_up, y_tr, z);
  stan::test::expect_ad(f_up, y_tr, y);
  stan::test::expect_ad_matvar(f, y, z);
  stan::test::expect_ad_matvar(f, u, y);
  stan::test::expect_ad_matvar(f_up, y_tr, z);
  stan::test::expect_ad_matvar(f_up, y_tr, y);

  Eigen::MatrixXd uu(2, 2);
  uu << 3, 0, 1, 4;
  Eigen::MatrixXd vv(2, 2);
  vv << 2, 3, 5, 7;
  stan::test::expect_ad(f, uu, vv);
  stan::test::expect_ad_matvar(f, uu, vv);

  // exception cases
  Eigen::MatrixXd d(3, 2);
  d << 1, 2, 5, -19, 73, 31;
  Eigen::VectorXd e(3);
  e << 1, 2, 9;
  stan::test::expect_ad(f, d, b);
  stan::test::expect_ad(f, d, c);
  stan::test::expect_ad(f, a, e);
  stan::test::expect_ad_matvar(f, d, b);
  stan::test::expect_ad_matvar(f, d, c);
  stan::test::expect_ad_matvar(f, a, e);

  Eigen::MatrixXd m33 = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd m44 = Eigen::MatrixXd::Zero(4, 4);
  Eigen::VectorXd v4 = Eigen::VectorXd::Zero(4);
  Eigen::RowVectorXd rv3 = Eigen::RowVectorXd::Zero(3);

  // exceptions: wrong sizes
  stan::test::expect_ad(f, m33, m44);
  stan::test::expect_ad(f, m33, v4);
  stan::test::expect_ad_matvar(f, m33, m44);
  stan::test::expect_ad_matvar(f, m33, v4);

  // exceptions: wrong types
  stan::test::expect_ad(f, m33, rv3);
  stan::test::expect_ad_matvar(f, m33, rv3);
}
