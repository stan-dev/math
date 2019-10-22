#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, divide) {
  auto f
      = [](const auto& x, const auto& y) { return stan::math::divide(x, y); };

  double x1 = 10;
  double x2 = -2;
  stan::test::expect_ad(f, x1, x2);

  Eigen::VectorXd v(1);
  v << 10;
  stan::test::expect_ad(f, v, x2);

  Eigen::RowVectorXd rv = v;
  stan::test::expect_ad(f, rv, x2);

  Eigen::MatrixXd m(1, 1);
  m << 10;
  stan::test::expect_ad(f, m, x2);

  Eigen::MatrixXd p(3, 2);
  p << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(f, p, x2);

  Eigen::VectorXd w(3);
  w << 100, 0, -3;
  stan::test::expect_ad(f, w, x2);

  Eigen::VectorXd u(4);
  u << 100, 0, -3, 4;
  stan::test::expect_ad(f, u, x2);
}
