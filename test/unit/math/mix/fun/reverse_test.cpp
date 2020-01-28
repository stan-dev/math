#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, reverse) {
  auto f = [](const auto& x) { return stan::math::reverse(x); };

  // 0 x 0
  Eigen::VectorXd x0(0);
  stan::test::expect_ad(f, x0);

  // 1 x 1
  Eigen::VectorXd x1(1);
  x1 << 1;
  stan::test::expect_ad(f, x1);

  // 4 x 4
  Eigen::VectorXd x(4);
  x << 1, 3, 4, -5;
  stan::test::expect_ad(f, x);
}
