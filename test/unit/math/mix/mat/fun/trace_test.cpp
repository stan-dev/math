#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, trace) {
  auto f = [](const auto& x) { return stan::math::trace(x); };

  Eigen::MatrixXd a(0, 0);
  stan::test::expect_ad(f, a);

  Eigen::MatrixXd b(1, 1);
  b << -1.2;
  stan::test::expect_ad(f, b);

  Eigen::MatrixXd c(2, 2);
  c << -1, 2, 5, 10;
  stan::test::expect_ad(f, a);

  Eigen::MatrixXd d(3, 3);
  d << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stan::test::expect_ad(f, d);
}
