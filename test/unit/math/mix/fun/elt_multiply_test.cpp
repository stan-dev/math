#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, elt_multiply_transpose_test) {
  auto f = [](const auto& x) {
    return stan::math::elt_multiply(x, x.transpose());
  };

  Eigen::MatrixXd x = Eigen::MatrixXd::Random(2, 2);

  stan::test::expect_ad_matvar(f, x);
}

TEST(MathMixMatFun, elt_multiply_block_test) {
  auto f = [](const auto& y, const auto& x) {
    return stan::math::elt_multiply(y.block(1, 1, 2, 2), x.block(0, 0, 2, 2));
  };

  Eigen::MatrixXd x(3, 3);

  stan::test::expect_ad_matvar(f, x, x);
}

TEST(MathMixMatFun, eltMultiply) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::elt_multiply(x, y);
  };

  double c = 1.1;
  Eigen::VectorXd a2(2);
  a2 << 2, 5;
  Eigen::VectorXd b2(2);
  b2 << 10, 100;
  stan::test::expect_ad(f, a2, b2);
  stan::test::expect_ad_matvar(f, a2, b2);
  stan::test::expect_ad(f, a2, c);
  stan::test::expect_ad_matvar(f, a2, c);
  stan::test::expect_ad(f, c, b2);
  stan::test::expect_ad_matvar(f, c, b2);

  Eigen::RowVectorXd c2(2);
  c2 << 2, 5;
  Eigen::RowVectorXd d2(2);
  d2 << 10, 100;
  stan::test::expect_ad(f, c2, d2);
  stan::test::expect_ad_matvar(f, c2, d2);
  stan::test::expect_ad(f, c2, c);
  stan::test::expect_ad_matvar(f, c2, c);
  stan::test::expect_ad(f, c, d2);
  stan::test::expect_ad_matvar(f, c, d2);

  Eigen::MatrixXd e23(2, 3);
  e23 << 2, 5, 7, 13, 29, 112;
  Eigen::MatrixXd f23(2, 3);
  f23 << -1.1, 2.2, 3.3, -1.4, 0.5, 6.0;
  stan::test::expect_ad(f, e23, f23);
  stan::test::expect_ad_matvar(f, e23, f23);
  stan::test::expect_ad(f, e23, c);
  stan::test::expect_ad_matvar(f, e23, c);
  stan::test::expect_ad(f, c, f23);
  stan::test::expect_ad_matvar(f, c, f23);

  Eigen::VectorXd v0(0);
  Eigen::RowVectorXd rv0(0);
  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, v0, v0);
  stan::test::expect_ad(f, rv0, rv0);
  stan::test::expect_ad(f, m00, m00);
  stan::test::expect_ad_matvar(f, v0, v0);
  stan::test::expect_ad_matvar(f, rv0, rv0);
  stan::test::expect_ad_matvar(f, m00, m00);
  stan::test::expect_ad(f, c, v0);
  stan::test::expect_ad(f, c, rv0);
  stan::test::expect_ad(f, c, m00);
  stan::test::expect_ad_matvar(f, c, v0);
  stan::test::expect_ad_matvar(f, c, rv0);
  stan::test::expect_ad_matvar(f, c, m00);
  stan::test::expect_ad(f, v0, c);
  stan::test::expect_ad(f, rv0, c);
  stan::test::expect_ad(f, m00, c);
  stan::test::expect_ad_matvar(f, v0, c);
  stan::test::expect_ad_matvar(f, rv0, c);
  stan::test::expect_ad_matvar(f, m00, c);
}
