#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, multiplyLowerTriSelfTranspose) {
  auto f = [](const auto& x) {
    return stan::math::multiply_lower_tri_self_transpose(x);
  };

  Eigen::MatrixXd a(0, 0);
  stan::test::expect_ad(f, a);

  Eigen::MatrixXd b(1, 1);
  b << 1.5;
  stan::test::expect_ad(f, b);

  Eigen::MatrixXd c(2, 2);
  c << 1, 0, 2, 3;
  stan::test::expect_ad(f, c);

  Eigen::MatrixXd d(3, 3);
  d << 1, 0, 0, 2, 3, 0, 4, 5, 6;
  stan::test::expect_ad(f, d);

  Eigen::MatrixXd e(3, 3);
  e << 1, 0, 1e6, 2, 3, 0, 4, 5, 6;
  stan::test::expect_ad(f, d);

  Eigen::MatrixXd g(2, 2);
  g << 3, 0, 4, -3;
  stan::test::expect_ad(f, g);

  Eigen::MatrixXd h = Eigen::MatrixXd::Zero(2, 3);
  stan::test::expect_ad(f, h);

  Eigen::MatrixXd j = Eigen::MatrixXd(2, 3);
  j << 1, 0, 0, 2, 3, 0;
  stan::test::expect_ad(f, j);
}
