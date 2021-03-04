#include <test/unit/math/test_ad.hpp>
#include <iostream>

void expect_diag_pre_multiply(const Eigen::VectorXd& v,
                              const Eigen::MatrixXd& a) {
  auto f = [](const auto& x, const auto& y) {
    return stan::math::diag_pre_multiply(x, y);
  };
  stan::test::expect_ad(f, v, a);
  stan::test::expect_ad_matvar(f, v, a);
  Eigen::RowVectorXd rv(v);
  stan::test::expect_ad(f, rv, a);
  stan::test::expect_ad_matvar(f, rv, a);
}
TEST(MathMixMatFun, diagPreMultiply) {
  using stan::test::relative_tolerance;
  // 0 x 0
  Eigen::MatrixXd a00(0, 0);
  Eigen::VectorXd u0(0);
  expect_diag_pre_multiply(u0, a00);

  // 1 x 1
  Eigen::MatrixXd a11(1, 1);
  a11 << 10;
  Eigen::VectorXd u1(1);
  u1 << 3;
  expect_diag_pre_multiply(u1, a11);

  // 2 x 2
  Eigen::MatrixXd a22(2, 2);
  a22 << 1, 10, 100, 1000;
  Eigen::VectorXd u2(2);
  u2 << 2, 3;
  expect_diag_pre_multiply(u2, a22);

  // 3 x 3
  Eigen::MatrixXd a33b(3, 3);
  a33b << 1, 2, 3, 2, 3, 4, 4, 5, 6;
  Eigen::VectorXd u3b(3);
  u3b << 1, 2, 3;
  expect_diag_pre_multiply(u3b, a33b);

  Eigen::MatrixXd a33c(3, 3);
  a33c << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Eigen::VectorXd u3c(3);
  u3c << 1, 2, 3;
  expect_diag_pre_multiply(u3c, a33c);

  Eigen::MatrixXd a33(3, 3);
  a33 << 1, 10, 100, 1000, 2, -4, 8, -16, 32;
  Eigen::VectorXd u3(3);
  u3 << -1.7, 111.2, -29.3;
  expect_diag_pre_multiply(u3, a33);

  Eigen::MatrixXd a33d(3, 3);
  a33d << 1, 0, 0, 0, 2, 0, 0, 0, 3;
  Eigen::VectorXd u3d(3);
  u3d << 1, 2, 3;
  expect_diag_pre_multiply(u3d, a33d);

  // error: mismatched sizes
  expect_diag_pre_multiply(u2, a33);
  expect_diag_pre_multiply(u3, a22);

  // non-square
  Eigen::MatrixXd b23(2, 3);
  b23 << 1, 2, 3, 4, 5, 6;
  expect_diag_pre_multiply(u2, b23);

  Eigen::MatrixXd b32(3, 2);
  b32 << 1, 2, 3, 4, 5, 6;
  expect_diag_pre_multiply(u3, b32);

  Eigen::MatrixXd b13(1, 3);
  b13 << 1, 2, 3;
  expect_diag_pre_multiply(u1, b13);

  Eigen::MatrixXd b31(3, 1);
  b31 << 1, 2, 3;
  expect_diag_pre_multiply(u3, b31);

  // non-square error: mismatched sizes
  expect_diag_pre_multiply(u3d, b23);
}
