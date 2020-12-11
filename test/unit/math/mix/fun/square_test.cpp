#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, square) {
  auto f = [](const auto& x1) { return stan::math::square(x1); };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -2.6, -1.0, -0.5, -0.2, 0.5, 1.3, 3, 5,
                                      1e5);

  Eigen::VectorXd x1(3);
  x1 << -1.0, 2.0, 3.0;
  Eigen::RowVectorXd x2(3);
  x2 << -1.0, 2.0, 3.0;
  Eigen::MatrixXd x3(2, 3);
  x3 << -1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
  stan::test::expect_ad_matvar(f, x1);
  stan::test::expect_ad_matvar(f, x2);
  stan::test::expect_ad_matvar(f, x3);

  Eigen::VectorXd x4(0);
  Eigen::RowVectorXd x5(0);
  Eigen::MatrixXd x6(0, 0);

  stan::test::expect_ad_matvar(f, x4);
  stan::test::expect_ad_matvar(f, x5);
  stan::test::expect_ad_matvar(f, x6);
}
