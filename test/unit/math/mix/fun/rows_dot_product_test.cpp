#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, rowsDotProduct) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::rows_dot_product(x, y);
  };

  Eigen::VectorXd v0(0);
  Eigen::RowVectorXd rv0(0);
  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, v0, v0);
  stan::test::expect_ad(f, rv0, rv0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad_matvar(f, v0, v0);
  stan::test::expect_ad_matvar(f, rv0, rv0);
  stan::test::expect_ad_matvar(f, m00, m00);

  Eigen::VectorXd u3(3);
  u3 << 1, 3, -5;
  Eigen::VectorXd v3(3);
  v3 << 4, -2, -1;
  stan::test::expect_ad(f, u3, v3);
  stan::test::expect_ad_matvar(f, u3, v3);

  Eigen::RowVectorXd ru3 = u3;
  Eigen::RowVectorXd rv3 = v3;
  stan::test::expect_ad(f, ru3, rv3);
  stan::test::expect_ad_matvar(f, ru3, rv3);

  Eigen::MatrixXd a33(3, 3);
  a33 << 1, 1, 1, 3, 3, 3, -5, -5, -5;
  Eigen::MatrixXd b33(3, 3);
  b33 << 4, 4, 4, -2, -2, -2, -1, -1, -1;
  stan::test::expect_ad(f, a33, b33);
  stan::test::expect_ad_matvar(f, a33, b33);

  Eigen::MatrixXd c32(3, 2);
  c32 << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd d32(3, 2);
  d32 << -1, -2, -3, -4, -5, -6;
  stan::test::expect_ad(f, c32, d32);
  stan::test::expect_ad_matvar(f, c32, d32);

  Eigen::MatrixXd c23 = c32.transpose();
  Eigen::MatrixXd d23 = d32.transpose();
  stan::test::expect_ad(f, c23, d23);
  stan::test::expect_ad_matvar(f, c23, d23);

  // exceptions---sizes
  Eigen::MatrixXd em33 = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd em32 = Eigen::MatrixXd::Zero(3, 2);
  Eigen::MatrixXd em23 = Eigen::MatrixXd::Zero(2, 3);

  Eigen::VectorXd ev3 = Eigen::VectorXd::Zero(3);
  Eigen::VectorXd ev2 = Eigen::VectorXd::Zero(2);

  Eigen::RowVectorXd erv3 = Eigen::RowVectorXd::Zero(3);
  Eigen::RowVectorXd erv2 = Eigen::RowVectorXd::Zero(2);

  stan::test::expect_ad(f, ev2, ev3);
  stan::test::expect_ad(f, erv2, erv3);
  stan::test::expect_ad(f, em33, em23);
  stan::test::expect_ad(f, em23, em33);
  stan::test::expect_ad_matvar(f, ev2, ev3);
  stan::test::expect_ad_matvar(f, erv2, erv3);
  stan::test::expect_ad_matvar(f, em33, em23);
  stan::test::expect_ad_matvar(f, em23, em33);
}
