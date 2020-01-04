#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, quadFormDiag) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::quad_form_diag(x, y);
  };
  Eigen::MatrixXd m00(0, 0);
  Eigen::VectorXd v0(0);

  Eigen::MatrixXd m11(1, 1);
  m11 << 1;
  Eigen::VectorXd v1(1);
  v1 << 2;

  Eigen::MatrixXd m22(2, 2);
  m22 << 2, 3, 4, 5;
  Eigen::VectorXd v2(2);
  v2 << 100, 10;

  Eigen::MatrixXd m33(3, 3);
  m33 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Eigen::VectorXd v3(3);
  v3 << 1, 2, 3;

  // matched sizes
  stan::test::expect_ad(f, m00, v0);
  stan::test::expect_ad(f, m11, v1);
  stan::test::expect_ad(f, m22, v2);
  stan::test::expect_ad(f, m33, v3);

  // exceptions from mismached sizes
  stan::test::expect_ad(f, m33, v2);
  stan::test::expect_ad(f, m22, v3);
}
