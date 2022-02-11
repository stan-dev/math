#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, symmetrize_from_lower_tri) {
  auto f
      = [](const auto& x) { return stan::math::symmetrize_from_lower_tri(x); };

  Eigen::MatrixXd a(0, 0);
  stan::test::expect_ad(f, a);

  Eigen::MatrixXd b(1, 1);
  b << -1.2;
  stan::test::expect_ad(f, b);

  Eigen::MatrixXd c(2, 2);
  c << -1, 2, 5, 10;
  stan::test::expect_ad(f, c);

  // Test invalid sizes
  Eigen::MatrixXd d(2, 3);
  d << -1, 2, -3, 5, 10, 100;
  stan::test::expect_ad(f, d);
}
