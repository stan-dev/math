#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, segment) {
  auto f = [](int i, int n) {
    return [=](const auto& y) { return stan::math::segment(y, i, n); };
  };

  Eigen::VectorXd a(0);
  stan::test::expect_ad(f(1, 0), a);  // out of bounds
  stan::test::expect_ad(f(1, 2), a);  // out of bounds

  Eigen::VectorXd b(1);
  b << 1.2;
  stan::test::expect_ad(f(1, 1), b);
  stan::test::expect_ad(f(1, 2), b);  // out of bounds

  Eigen::VectorXd v(3);
  v << 1, 2, 3;
  Eigen::RowVectorXd rv = v;
  stan::test::expect_ad(f(1, 0), v);
  stan::test::expect_ad(f(1, 3), v);
  stan::test::expect_ad(f(1, 0), rv);
  stan::test::expect_ad(f(1, 3), rv);

  // exception (out of bounds)
  stan::test::expect_ad(f(1, 4), v);   // out of bounds
  stan::test::expect_ad(f(1, 4), rv);  // out of bounds
  stan::test::expect_ad(f(1, 7), v);   // out of bounds
  stan::test::expect_ad(f(1, 7), rv);  // out of bounds

  Eigen::VectorXd u(4);
  u << 1, 2, 3, 4;
  Eigen::RowVectorXd ru = u;
  stan::test::expect_ad(f(2, 2), ru);
}
