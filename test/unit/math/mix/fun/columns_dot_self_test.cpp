#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, columnsDotSelf) {
  auto f = [](const auto& x) { return stan::math::columns_dot_self(x); };
  Eigen::MatrixXd a11(1, 1);
  a11 << 2;
  stan::test::expect_ad(f, a11);
  stan::test::expect_ad_matvar(f, a11);

  Eigen::MatrixXd a12(1, 2);
  a12 << 2, 3;
  stan::test::expect_ad(f, a12);
  stan::test::expect_ad_matvar(f, a12);

  Eigen::MatrixXd a22(2, 2);
  a22 << 2, 3, 4, 5;
  stan::test::expect_ad(f, a22);
  stan::test::expect_ad_matvar(f, a22);

  Eigen::VectorXd u3(3);
  u3 << 1, 3, -5;
  Eigen::VectorXd v3(3);
  v3 << 4, -2, -1;
  stan::test::expect_ad(f, u3);
  stan::test::expect_ad_matvar(f, u3);

  Eigen::RowVectorXd ru3 = u3;
  stan::test::expect_ad(f, ru3);
  stan::test::expect_ad_matvar(f, ru3);

  Eigen::MatrixXd a33(3, 3);
  a33 << 1, 1, 1, 3, 3, 3, -5, -5, -5;
  Eigen::MatrixXd b33(3, 3);
  b33 << 4, 4, 4, -2, -2, -2, -1, -1, -1;
  stan::test::expect_ad(f, a33);
  stan::test::expect_ad(f, b33);
  stan::test::expect_ad_matvar(f, a33);
  stan::test::expect_ad_matvar(f, b33);

  Eigen::MatrixXd c32(3, 2);
  c32 << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd d32(3, 2);
  d32 << -1, -2, -3, -4, -5, -6;
  stan::test::expect_ad(f, c32);
  stan::test::expect_ad(f, d32);
  stan::test::expect_ad_matvar(f, c32);
  stan::test::expect_ad_matvar(f, d32);

  Eigen::MatrixXd c23 = c32.transpose();
  Eigen::MatrixXd d23 = d32.transpose();
  stan::test::expect_ad(f, c23);
  stan::test::expect_ad(f, d23);
  stan::test::expect_ad_matvar(f, c23);
  stan::test::expect_ad_matvar(f, d23);

  Eigen::MatrixXd a00(0, 0);
  Eigen::VectorXd v0(0);
  Eigen::RowVectorXd rv0(0);
  stan::test::expect_ad(f, a00);
  stan::test::expect_ad(f, v0);
  stan::test::expect_ad(f, rv0);
  stan::test::expect_ad_matvar(f, a00);
  stan::test::expect_ad_matvar(f, v0);
  stan::test::expect_ad_matvar(f, rv0);
}
