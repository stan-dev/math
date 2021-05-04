#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, distance) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::distance(x, y); };

  // 0 x 0
  Eigen::VectorXd x0(0);
  Eigen::VectorXd y0(0);
  stan::test::expect_ad(f, x0, y0);
  stan::test::expect_ad_matvar(f, x0, y0);

  // 1 x 1
  Eigen::VectorXd x1(1);
  x1 << 1;
  Eigen::VectorXd y1(1);
  y1 << -2.3;
  stan::test::expect_ad(f, x1, y1);
  stan::test::expect_ad_matvar(f, x1, y1);

  // 2 x 2
  Eigen::VectorXd x2(2);
  x2 << 2, -3;
  Eigen::VectorXd y2(2);
  y2 << -2.3, 1.1;
  stan::test::expect_ad(f, x2, y2);
  stan::test::expect_ad_matvar(f, x2, y2);

  // 3 x 3
  Eigen::VectorXd x(3);
  x << 1, 3, -5;
  Eigen::VectorXd y(3);
  y << 4, -2, -1;
  stan::test::expect_ad(f, x, y);
  stan::test::expect_ad_matvar(f, x, y);

  // exception cases
  Eigen::VectorXd z(2);
  z << 1, 2;
  stan::test::expect_ad(f, x, z);
  stan::test::expect_ad(f, z, x);
  stan::test::expect_ad_matvar(f, x, z);
  stan::test::expect_ad_matvar(f, z, x);
}
