#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, trace) {
  auto f = [](const auto& x) { return stan::math::trace(x); };

  Eigen::MatrixXd a(0, 0);
  stan::test::expect_ad(f, a);
  stan::test::expect_ad_matvar(f, a);

  Eigen::MatrixXd b(1, 1);
  b << -1.2;
  stan::test::expect_ad(f, b);
  stan::test::expect_ad_matvar(f, b);

  Eigen::MatrixXd c(2, 2);
  c << -1, 2, 5, 10;
  stan::test::expect_ad(f, c);
  stan::test::expect_ad_matvar(f, c);

  Eigen::MatrixXd d(3, 3);
  d << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::test::expect_ad(f, d);
  stan::test::expect_ad_matvar(f, d);

  Eigen::MatrixXd e(2, 3);
  e << -1, 2, 5, 10, 1.0, 1.5;
  stan::test::expect_ad(f, e);
  stan::test::expect_ad_matvar(f, e);

  Eigen::MatrixXd h(3, 2);
  h << -1, 2, 5, 10, 1.0, 1.5;
  stan::test::expect_ad(f, h);
  stan::test::expect_ad_matvar(f, h);
}
