#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, block) {
  auto f = [](int i, int j, int m, int n) {
    return [=](const auto& y) { return stan::math::block(y, i, j, m, n); };
  };
  Eigen::MatrixXd b(1, 1);
  b << 3.2;
  stan::test::expect_ad(f(1, 1, 1, 1), b);

  Eigen::MatrixXd a(3, 3);
  a << 1, 4, 9, 1, 4, 9, 1, 4, 9;
  stan::test::expect_ad(f(1, 1, 3, 3), a);
  stan::test::expect_ad(f(2, 2, 2, 2), a);

  // exceptions (out of range)
  stan::test::expect_ad(f(0, 0, 1, 1), a);
  stan::test::expect_ad(f(0, 0, 4, 4), a);
}
