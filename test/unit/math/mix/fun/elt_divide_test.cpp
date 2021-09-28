#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, elt_divide_transpose_test) {
  auto f
      = [](const auto& x) { return stan::math::elt_divide(x, x.transpose()); };

  Eigen::MatrixXd x(2, 2);

  stan::test::expect_ad_matvar(f, x);
}

TEST(MathMixMatFun, elt_divide_block_test) {
  auto f = [](const auto& y, const auto& x) {
    return stan::math::elt_divide(y.block(1, 1, 2, 2), x.block(0, 0, 2, 2));
  };

  Eigen::MatrixXd x(3, 3);

  stan::test::expect_ad_matvar(f, x, x);
}

TEST(MathMixMatFun, eltDivide) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::elt_divide(x, y);
  };

  Eigen::VectorXd a2(2);
  a2 << 2, 5;
  Eigen::VectorXd b2(2);
  b2 << 10, 100;
  stan::test::expect_ad(f, a2, b2);
  stan::test::expect_ad_matvar(f, a2, b2);

  Eigen::RowVectorXd c2(2);
  c2 << 2, 5;
  Eigen::RowVectorXd d2(2);
  d2 << 10, 100;
  stan::test::expect_ad(f, c2, d2);
  stan::test::expect_ad_matvar(f, c2, d2);

  Eigen::MatrixXd e23(2, 3);
  e23 << 2, 5, 7, 13, 29, 112;
  Eigen::MatrixXd f23(2, 3);
  f23 << 1e1, 1e2, 1e3, 1e4, 1e5, 1e6;
  stan::test::expect_ad(f, e23, f23);
  stan::test::expect_ad_matvar(f, e23, f23);

  double z = 10;
  stan::test::expect_ad(f, a2, z);
  stan::test::expect_ad(f, c2, z);
  stan::test::expect_ad(f, e23, z);

  stan::test::expect_ad_matvar(f, a2, z);
  stan::test::expect_ad_matvar(f, c2, z);
  stan::test::expect_ad_matvar(f, e23, z);

  stan::test::expect_ad(f, z, a2);
  stan::test::expect_ad(f, z, c2);
  stan::test::expect_ad(f, z, e23);

  stan::test::expect_ad_matvar(f, z, a2);
  stan::test::expect_ad_matvar(f, z, c2);
  stan::test::expect_ad_matvar(f, z, e23);

  Eigen::VectorXd v0(0);
  Eigen::RowVectorXd rv0(0);
  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, v0, v0);
  stan::test::expect_ad(f, rv0, rv0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad(f, v0, z);
  stan::test::expect_ad(f, rv0, z);
  stan::test::expect_ad(f, m00, z);

  stan::test::expect_ad_matvar(f, v0, v0);
  stan::test::expect_ad_matvar(f, rv0, rv0);
  stan::test::expect_ad_matvar(f, m00, m00);
  stan::test::expect_ad_matvar(f, v0, z);
  stan::test::expect_ad_matvar(f, rv0, z);
  stan::test::expect_ad_matvar(f, m00, z);
}
