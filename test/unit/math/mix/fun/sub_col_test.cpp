#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, subCol) {
  auto f = [](int i, int j, int k) {
    return [=](const auto& y) { return stan::math::sub_col(y, i, j, k); };
  };

  Eigen::MatrixXd a(1, 1);
  a << 3.2;
  stan::test::expect_ad(f(1, 1, 0), a);
  stan::test::expect_ad(f(1, 1, 1), a);

  Eigen::MatrixXd b(3, 4);
  b << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;
  stan::test::expect_ad(f(1, 1, 0), b);
  stan::test::expect_ad(f(1, 1, 1), b);
  stan::test::expect_ad(f(1, 1, 3), b);
  stan::test::expect_ad(f(1, 2, 2), b);
  stan::test::expect_ad(f(3, 4, 1), b);
  stan::test::expect_ad(f(2, 3, 2), b);
  stan::test::expect_ad(f(1, 1, 7), b);  // exception--range
  stan::test::expect_ad(f(7, 1, 1), b);  // exception--range
  stan::test::expect_ad(f(1, 7, 1), b);  // exception--range

  Eigen::MatrixXd c(0, 0);
  stan::test::expect_ad(f(0, 0, 0), c);
  stan::test::expect_ad(f(0, 1, 0), c);
  stan::test::expect_ad(f(0, 1, 1), c);
  stan::test::expect_ad(f(1, 0, 0), c);
  stan::test::expect_ad(f(1, 0, 1), c);
  stan::test::expect_ad(f(1, 1, 0), c);
  stan::test::expect_ad(f(1, 1, 1), c);
}
